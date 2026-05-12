# Source code for matlantis_features.features.md.md
    
    
    import logging
    import os
    import tempfile
    from dataclasses import dataclass
    from pathlib import Path
    from shutil import copyfile
    from subprocess import check_output
    from typing import Any, Dict, List, Optional, Tuple
    
    import numpy as np
    from ase import units
    from ase.io.trajectory import Trajectory, TrajectoryReader
    from tqdm.auto import tqdm
    from tqdm.contrib.logging import logging_redirect_tqdm
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    from .ase_integrators import MDIntegratorBase
    from .ase_simulation import ASECheckpoint, ASESimulation, ASETrajWriter
    from .md_extensions import MDExtensionBase
    from .md_simulation_base import MDExtensionList
    from .md_system_base import MDSystemBase
    
    logger = logging.getLogger(__name__)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult)@dataclass
    class MDFeatureResult:
        """A dataclass for result of MDFeature."""
    
        traj_path: Optional[str] = None
        checkpoint_path: Optional[str] = None
        temp_dir: Optional[tempfile.TemporaryDirectory] = None  # type: ignore
    
        # TODO: remove this compatibility method
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult.get_traj_obj)    def get_traj_obj(self, ctxt_save_dir: Path) -> TrajectoryReader:
            """Get ASE's trajectory reader object of the MD result.
    
            Args:
                ctxt_save_dir (Path):
                  Context's temporary directory location.
                  Ignored if this obj contains trajectory path in absolute path.
            Returns:
                TrajectoryReader :
                  TrajectoryReader object for this MD result.
            """
            if self.traj_path is None:
                raise RuntimeError("traj is not saved for this md result")
            path = Path(self.traj_path)
            if path.is_absolute():
                return Trajectory(path, "r")
            else:
                return Trajectory(ctxt_save_dir / path, "r")
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult.save)    def save(self, traj_file_name: Optional[str] = None, checkpoint_file_name: Optional[str] = None) -> None:
            """Save the MD trajectory and checkpoint files to a different location.
    
            Args:
                traj_file_name (str or None, optional):
                  The new location of MD trajectory file. Defaults to None.
                checkpoint_file_name (str or None, optional):
                  The new location of MD checkpoint file. Defaults to None.
            """
            if traj_file_name:
                result_fname = traj_file_name if Path(traj_file_name).is_absolute() else str(Path.cwd() / traj_file_name)
                copyfile(self.traj_path, result_fname)  # type: ignore
            if checkpoint_file_name:
                result_fname = checkpoint_file_name if Path(checkpoint_file_name).is_absolute() else str(Path.cwd() / checkpoint_file_name)
                copyfile(self.checkpoint_path, result_fname)  # type: ignore
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult.from_traj_obj)    @classmethod
        def from_traj_obj(cls, traj_file_name: str) -> "MDFeatureResult":
            """Get MDFeatureResult from trajectory file.
    
            Args:
                traj_file_name (str): Path to the trajectory file.
            Returns:
                MDFeatureResult : MDFeatureResult.
            """
            path = Path(traj_file_name).absolute()
            if not path.is_file():
                raise FileNotFoundError
            traj = Trajectory(str(path), "r")
            if "time" not in traj[0].info or "step" not in traj[0].info:
                logger.warning("Time information is missing, please make sure the trajectory is generated from MDFeature.")
            return cls(traj_path=str(path))
    
    
    
    
    class _ProgressExtension(MDExtensionBase):
        """Extension to display progress bar and log information."""
    
        def __init__(
            self,
            md_feature: "MDFeature",
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 100,
        ):
            """Create the extension.
    
            Args:
                md_feature (MDFeature): Target MD feature object.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 10.
            """
            self.md_feature = md_feature
            ntotal = self.md_feature.n_run
            self.show_progress_bar = show_progress_bar
            if show_progress_bar:
                if tqdm_options is None:
                    self._pbar = tqdm(total=ntotal)
                else:
                    self._pbar = tqdm(total=ntotal, **tqdm_options)
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            if show_logger:
                self._logger = logging.getLogger(__name__)
    
        def __call__(self, system: MDSystemBase, integrator: MDIntegratorBase) -> None:
            """Call the feature's report method.
    
            Args:
                system (MDSystemBase):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            # pbar is not updated at the init frame (step==0)
            with logging_redirect_tqdm():
                if self.show_logger:
                    if system.current_step % self.logger_interval == 0:
                        atoms = system.ase_atoms  # type: ignore
                        if not np.all(atoms.get_pbc()):
                            volume = 0.0
                            density = 0.0
                        else:
                            volume = atoms.get_volume()
                            density = atoms.get_masses().sum() / atoms.get_volume() / units.kg * 1e27
                        logger.info(
                            f"steps: {system.current_step:5d}  "
                            f"energy：{atoms.get_potential_energy() / len(atoms):5.2f} eV/atom  "
                            f"total energy: {atoms.get_total_energy() / len(atoms):5.2f} eV/atom  "
                            f"temperature: {atoms.get_temperature():5.2f} K  "
                            f"volume: {volume:5.0f} Ang^3  "
                            f"density: {density:5.3f} g/cm^3"
                        )
                if self.show_progress_bar and system.current_step > 0:
                    self._pbar.update(1)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeature.html#matlantis_features.features.md.md.MDFeature)class MDFeature(FeatureBase):
        """The matlantis-feature for molecular dynamics simulation."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeature.html#matlantis_features.features.md.md.MDFeature.__init__)    def __init__(
            self,
            integrator: MDIntegratorBase,
            n_run: int,
            checkpoint_file_name: Optional[str] = None,
            checkpoint_freq: Optional[int] = 10,
            traj_file_name: Optional[str] = None,
            traj_freq: Optional[int] = 1,
            traj_props: Optional[List[str]] = None,
            traj_append: bool = False,
            traj_init_frame: bool = False,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 100,
            estimator_fn: Optional[EstimatorFnType] = None,
        ):
            """Initialize an instance.
    
            Args:
                integrator (MDIntegratorBase): Integrator object used for the simulation.
                n_run (int): Number or the steps to be run.
                checkpoint_file_name (str or None, optional):
                  Checkpoint file name. If this value is None and checkpoint_freq>0,
                  checkpoint file name is automatically generated.
                  If relative path name is given,
                  the file is stored in the current working directory.
                  Defaults to None.
                checkpoint_freq (int or None, optional):
                  Checkpoint creation frequency.
                  If this value is None or non positive, no checkpoit file is created.
                  Defaults to 10.
                traj_file_name (str or None, optional):
                  Trajectory file name. If this value is None and traj_freq>0,
                  trajectory file name is automatically generated.
                  If relative path name is given,
                  the file is stored in the current working directory.
                  Defaults to None.
                traj_freq (int or None, optional):
                  Trajectory saving frequency.
                  If this value is None or non positive, no trajectory file is created.
                  Defaults to 1.
                traj_props (list[str] or None, optional):
                  Names of the properties stored in the trajectory file.
                  If None, energy, forces, and stress properties are stored.
                  Defaults to None.
                traj_append (bool, optional):
                  If True and trajectory file exists, trajectory will be
                  appended to the existing file. Defaults to False.
                traj_init_frame (bool, optional):
                  If True, the initial structure of the simulation is
                  also stored in the trajectory file. Defaults to False.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 100.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            super().__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            with self.init_scope():
                self.integrator = integrator
                self.n_run = n_run
    
                self.checkpoint_file = checkpoint_file_name
                self.checkpoint_freq = checkpoint_freq
    
                self.trajectory_file = traj_file_name
                self.traj_freq = traj_freq
                if traj_props is None:
                    self.traj_props = ["energy", "forces", "stress"]
                else:
                    self.traj_props = traj_props
                if traj_append:
                    self.traj_write_mode = "a"
                else:
                    self.traj_write_mode = "w"
                self.traj_init_frame = traj_init_frame
                self._temp_traj_dir: Optional[tempfile.TemporaryDirectory] = None  # type: ignore
                self.show_pbar = show_progress_bar
                self.tqdm_options = tqdm_options
                self.show_logger = show_logger
                self.logger_interval = logger_interval
                self.estimator_fn = estimator_fn
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md.MDFeature.html#matlantis_features.features.md.md.MDFeature.get_temp_dir)    def get_temp_dir(self) -> Path:
            """Create temperory directory.
    
            Returns:
                Path : The location of temperory directory.
            """
            if self._temp_traj_dir is None:
                self._temp_traj_dir = tempfile.TemporaryDirectory(prefix="matlantis_")
            return Path(self._temp_traj_dir.name)
    
    
    
        def _check_disk_size(self, path: str) -> Tuple[str, str]:
            dir = Path(path).parent
            info = check_output(["df", "-h", dir]).split()
            return info[-1].decode("UTF-8"), info[-5].decode("UTF-8")
    
        def _make_filename(self, fname: Optional[str], suffix: str) -> str:
            result_fname = None
            if fname is None:
                handle, result_fname = tempfile.mkstemp(suffix=suffix, dir=self.get_temp_dir())
                os.close(handle)
            else:
                if Path(fname).is_absolute():
                    result_fname = fname
                else:
                    result_fname = str(Path.cwd() / fname)
            if not os.path.isdir(os.path.dirname(result_fname)):
                os.makedirs(os.path.dirname(result_fname), exist_ok=True)
            if suffix == ".traj":
                mount, size = self._check_disk_size(result_fname)
                logger.warning(f"The MD trajectory will be saved at {result_fname}.")
                logger.warning(f"Note: The max disk size of {mount} is about {size}.")
                if fname is None:
                    logger.warning("Note: The MD trajectory is saved in a temporary directory, which will be automatically deleted later.")
            return result_fname
    
        def __call__(self, system: MDSystemBase, extensions: Optional[MDExtensionList] = None) -> MDFeatureResult:
            """Run the MD simulation.
    
            Args:
                system (MDSystemBase):
                  Target system of the MD simulation.
                extensions (MDExtensionList or None, optional):
                  Extension list (tuple of extension object and interval) used for the simulation.
                  Defaults to None.
            Returns:
                MDFeatureResult : Result object for the MD simulation.
            """
            sim = ASESimulation(system, self.integrator, estimator_fn=self.estimator_fn)
    
            if self.show_pbar or self.show_logger:
                sim.append_extension(
                    _ProgressExtension(self, self.show_pbar, self.tqdm_options, self.show_logger, self.logger_interval),
                    1,
                )
    
            # Add custom extensions
            #   Custom extensions should be added before the checkpoint/trajectory exts,
            #   because the checkpoint/trajectory exts have to save the information
            #   calculated by the custom extensions.
            if extensions is not None:
                for extn in extensions:
                    sim.append_extension(*extn)
    
            # Add checkpoint extension
            checkpoint_fname = None
            if self.checkpoint_freq is not None and self.checkpoint_freq > 0:
                checkpoint_fname = self._make_filename(self.checkpoint_file, suffix=".chkp")
                chkpw = ASECheckpoint(file_name=checkpoint_fname)
                sim.append_extension(chkpw, self.checkpoint_freq)
    
            # Add trajectory extension
            trajw = None
            traj_fname = None
            if self.traj_freq is not None and self.traj_freq > 0:
                traj_fname = self._make_filename(self.trajectory_file, suffix=".traj")
                trajw = ASETrajWriter(
                    file_name=traj_fname,
                    write_mode=self.traj_write_mode,
                    props=self.traj_props,
                    save_init_frame=self.traj_init_frame,
                )
                sim.append_extension(trajw, self.traj_freq)
    
            sim.run(n_run_steps=self.n_run)
    
            if trajw is not None:
                trajw.close()
    
            return MDFeatureResult(
                # md_feature=self,
                traj_path=traj_fname,
                checkpoint_path=checkpoint_fname,
                temp_dir=self._temp_traj_dir,
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=self.n_run, time_additional=0.1)
    
    
    
