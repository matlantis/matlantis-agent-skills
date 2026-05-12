# Source code for matlantis_features.features.reaction.rest_scan
    
    
    import logging
    import pathlib
    from dataclasses import dataclass
    from logging import Logger
    from typing import Any, Callable, Dict, List, Optional, Union
    
    import numpy as np
    import plotly.graph_objects as go
    from ase import Atoms
    from ase.io import Trajectory
    from restscan import (
        ExitByEnergy,
        ExitedByEnergy,
        Scan,
        ScanDihedral,  # noqa: F401 ; not used by imported for docs
        ScanDistance,  # noqa: F401 ; not used by imported for docs
        ScanRestraint,
        TrajectoryLogger,
        calculate_dihedral,  # noqa: F401 ; not used by imported for docs
    )
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.utils import get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult)@dataclass
    class RestScanFeatureResult:
        """A dataclass for result of RestScanFeature."""
    
        trajectory_logger: TrajectoryLogger
        traj_path: Optional[str]
        converged: bool
        warnflag: Optional[int]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.from_trajectory)    @classmethod
        def from_trajectory(cls, traj_file_name: Union[str, pathlib.Path]) -> "RestScanFeatureResult":
            """Get RestScanFeatureResult from trajectory file.
    
            Args:
                traj_file_name (str or pathlib.Path): Path to the trajectory file.
            Returns:
                RestScanFeatureResult : RestScanFeatureResult.
            """
            if isinstance(traj_file_name, str):
                traj_path = str(pathlib.Path(traj_file_name).absolute())
            elif isinstance(traj_file_name, pathlib.Path):
                traj_path = str(traj_file_name.absolute())
            else:
                raise ValueError
    
            traj = [atoms for atoms in Trajectory(traj_path)]
            images = TrajectoryLogger.summarize_static(traj, kink_dx=0.0, local_minima_de=None)
            trajectory_logger = TrajectoryLogger(atoms=images[-1])
            trajectory_logger.images = images
    
            return cls(trajectory_logger=trajectory_logger, traj_path=traj_path, converged=True, warnflag=None)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.get_images)    def get_images(self, kink_dx: float = 0.1, local_minima_de: Optional[float] = None) -> List[Atoms]:
            """Summerize the images of the restscan result.
    
            Args:
                kink_dx (float, optional):
                    Images whose maximum distance is not more than dx are regarded as identical
                    and omitted. Defaults to 0.1.
                local_minima_de (float or None, optional):
                    If de is given, returns a summary up to the last local minimum.
                    Here, de is the threshold to be regarded as local minimum. Defaults to None.
            Returns:
                List[Atoms]: A list of images.
            """
            images: List[Atoms] = self.trajectory_logger.summarize(kink_dx=kink_dx, local_minima_de=local_minima_de)
            return images
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.plot)    def plot(
            self,
            plt_name: Optional[str] = None,
            kink_dx: float = 0.1,
            local_minima_de: Optional[float] = None,
        ) -> go.Figure:
            """Plot the energy profile along the images of the restscan result.
    
            Args:
                plt_name (str or None, optional):
                    The file name to which the plot of the energy profile
                    will be saved. Supported formats include 'png' 'jpg' 'jpeg' 'webp' 'svg'
                    and 'pdf'. If 'None' is provided, the figure will not be saved. Defaults to None.
                kink_dx (float, optional):
                    Images whose maximum distance is not more than dx are regarded as identical
                    and omitted. Defaults to 0.1.
                local_minima_de (float or None, optional):
                    If de is given, returns a summary up to the last local minimum.
                    Here, de is the threshold to be regarded as local minimum. Defaults to None.
            Returns:
                go.Figure:
                    The plot of the energy profile in the restscan calculation result as a
                    plotly.graph_objects.Figure instance.
            """
            images = self.get_images(kink_dx=kink_dx, local_minima_de=local_minima_de)
            energies = [atoms.get_potential_energy() for atoms in images]
            steps = np.arange(len(energies))
            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=energies,
                    showlegend=False,
                )
            )
            fig.update_xaxes(title="reaction coordination")
            fig.update_yaxes(title="energy (eV)")
            if plt_name is not None:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeature.html#matlantis_features.features.reaction.rest_scan.RestScanFeature)class RestScanFeature(FeatureBase):
        """RestScan is a method for continuously changing molecular structures into different structures.\
        This can be used to make initial estimates of reaction paths or to obtain alternative structures \
        from existing molecular structures. It's similar to the Gaussian scan function.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/rest_scan.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.rest_scan.RestScanFeature.html#matlantis_features.features.reaction.rest_scan.RestScanFeature.__init__)    def __init__(
            self,
            scans: List[List[ScanRestraint]],
            trajectory: Optional[str] = None,
            logfile: Optional[str] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            estimator_fn: Optional[EstimatorFnType] = None,
        ):
            """Initialize an instance.
    
            Args:
                scans (list[list[ScanRestraint]]):
                    The scans is given as a list of lists of scans.
                    The number of scans to be applied is increased in the order of the outer list.
                trajectory (str or None, optional):
                    Attach trajectory object.  If *trajectory* is a string a
                    Trajectory will be constructed.  Use *None* for no
                    trajectory. Defaults to None.
                logfile (str or None, optional):
                    If logfile is a string, a file with that name will be opened.
                    Use '-' for stdout. Defaults to None.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 10.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            logger: Optional[Logger]
            if show_logger:
                logger = logging.getLogger(__name__)
            else:
                logger = None
    
            def get_scan(atoms: ASEAtoms) -> Scan:
                scan = Scan(
                    atoms=atoms,
                    scans=scans,
                    logfile=logfile,
                    trajectory=trajectory,
                    show_progress_bar=show_progress_bar,
                    tqdm_options=tqdm_options,
                    logger=logger,
                    logger_interval=logger_interval,
                )
                return scan
    
            super(RestScanFeature, self).__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            with self.init_scope():
                self.get_scan = get_scan
                self.scans = scans
                self.trajectory = trajectory
                self.logfile = logfile
                self.show_progress_bar = show_progress_bar
                self.tqdm_option = tqdm_options
                self.show_logger = show_logger
                self.logger_interval = logger_interval
                self.estimator_fn = estimator_fn
    
    
    
        def __call__(
            self,
            atoms: Union[ASEAtoms, MatlantisAtoms],
            fmax: float = 0.01,
            n_run: int = 10000,
            energy_stop: Optional[float] = None,
            local_minima_de: float = 0.01,
            local_minima_n: int = 1,
            attach_methods: Optional[List[Callable[[], None]]] = None,
        ) -> RestScanFeatureResult:
            """Run the RestScan calculation.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms): The input structure.
                fmax (float, optional):
                    The maximum force (in eV/Angstrom) for convergence of scan.
                    Defaults to 0.01.
                n_run (int, optional):
                    The maximum number of scan steps performed to try to
                    reach force convergence (specified by fmax). Defaults to 10000.
                energy_stop (float or None, optional):
                    Exit when the potential energy of the atoms is higher than this value.
                    Defaults to None.
                local_minima_de (float, optional):
                    Local minimas that are lower than local maxima by this value
                    are counted as local minimas. Defaults to 0.01.
                local_minima_n (int, optional):
                    Exit when the trajectory reach n'th local minima. Defaults to 1.
                attach_methods (list[ -> None] or None, optional):
                    Functions to be called by the restscan. A list of callable functions.
                    If None, no functions will be attached. Defaults to None.
            Returns:
                RestScanFeatureResult : The calculation result.
            """
            system: ASEAtoms
            if isinstance(atoms, MatlantisAtoms):
                system = atoms.ase_atoms
            else:
                system = atoms
    
            scan = self.get_scan(system)
    
            with get_calculator(estimator_fn=self.estimator_fn) as calculator:
                system.calc = calculator
    
                trajectory_logger = TrajectoryLogger(system)
                scan.attach(trajectory_logger)
    
                if energy_stop is not None:
                    exit = ExitByEnergy(
                        system,
                        energy_stop=energy_stop,
                        local_minima_de=local_minima_de,
                        local_minima_n=local_minima_n,
                    )
                    scan.attach(exit)
    
                if attach_methods is not None:
                    for method in attach_methods:
                        scan.attach(method)
    
                warnflag: Optional[int] = None
                try:
                    scan.run(fmax=fmax, steps=n_run)
                except ExitedByEnergy as e:
                    if self.show_logger:
                        logger = logging.getLogger(__name__)
                        logger.warning(e.args[0])
                    if "local minima reached" in e.args[0]:
                        converged = True
                    elif "higher than energy_stop" in e.args[0]:
                        converged = True
                        warnflag = 2
                    else:
                        raise NotImplementedError
                else:
                    converged = scan.converged()
                    if not converged:
                        warnflag = 1
    
            trajectory: Optional[str]
            if self.trajectory is not None:
                trajectory = str(pathlib.Path(self.trajectory).absolute())
            else:
                trajectory = self.trajectory
    
            return RestScanFeatureResult(
                trajectory_logger=trajectory_logger,
                traj_path=trajectory,
                converged=converged,
                warnflag=warnflag,
            )
    
    
    
