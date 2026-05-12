# Source code for matlantis_features.features.reaction.reaction_string
    
    
    import warnings
    from copy import deepcopy
    from dataclasses import dataclass
    from datetime import timedelta
    from pathlib import Path
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union
    
    import ase
    import matplotlib.pyplot as plt
    import numpy as np
    import plotly.graph_objects as go
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from ase.io.trajectory import Trajectory
    from reactionstring import DumpImages, ReactionPath
    from reactionstring.optimizers import FIRELBFGS, FIRES
    from reactionstring.utils.typing import OptimizerBuilder
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase, FeatureBaseOutput, FeatureBaseResult
    from matlantis_features.utils import get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    from matlantis_features.utils.save_util import (
        MatlantisNonSerializable,
        get_calculator_info,
        get_user_info,
        get_version_info,
        to_serialized_dict,
    )
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput)@dataclass
    class ReactionStringFeatureOutput(FeatureBaseOutput):
        """A dataclass containing the output of ReactionStringFeature."""
    
        num_segments: int
        reaction_string_images: List[List[MatlantisAtoms]]
        transition_state_indices: List[Optional[int]]
        imaginary_eigenvectors: List[np.ndarray]
        converged: bool
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult)@dataclass
    class ReactionStringFeatureResult(FeatureBaseResult):
        """A dataclass containing the result for ReactionStringFeature."""
    
        feature: "ReactionStringFeature"
        output: "ReactionStringFeatureOutput"
    
        def _flatten_images(self) -> Tuple[List[MatlantisAtoms], List[int], List[int]]:
            # Concatenate all images into a flat array
            #
            # Input:
            # [l1, ..., r1], [l2, ..., r2], ..., [ln, rn]
            # where (r1 = l2)
            #
            # Output:
            # [l1, ..., {r1 - 1}, l2, ... {r2 - 1}, ..., rn]
            # s.t. removing duplicate elements
            offset = 0
            flatten_atoms: List[MatlantisAtoms] = []
            stable_point_idx: List[int] = []
            saddle_point_idx: List[int] = []
    
            for images, transition_state_idx_in_segment in zip(
                self.output.reaction_string_images,
                self.output.transition_state_indices,
            ):
                flatten_atoms.extend(images[:-1])
                stable_point_idx.append(offset)
                if transition_state_idx_in_segment is not None:
                    saddle_point_idx.append(offset + transition_state_idx_in_segment)
                offset += len(images) - 1
    
            # Above loop tracks left stable points of each peeks.
            # For the last peek, we need to add the right stable point manually.
            flatten_atoms.append(self.output.reaction_string_images[-1][-1])
            stable_point_idx.append(offset)
            return flatten_atoms, stable_point_idx, saddle_point_idx
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.plot_energy)    def plot_energy(self, segment_id: Optional[int] = None, plt_name: Optional[str] = None) -> go.Figure:
            """Plot the energy profile along the reaction path.
    
            Args:
                segment_id (int or None, optional):
                    The index of the segment of path to plot. If None, all the path segments will be
                    plotted togather. Defaults to None.
                plt_name (str or None, optional):
                    The file name to which the plot of the energy profile
                    will be saved. Supported formats include 'png' 'jpg' 'jpeg' 'webp' 'svg'
                    and 'pdf'. If 'None' is provided, the figure will not be saved. Defaults to None.
            Returns:
                go.Figure :
                    The plot of the energy profile in the reaction string calculation result as a
                    plotly.graph_objects.Figure instance.
            """
            if segment_id is not None:
                flatten_atoms = self.output.reaction_string_images[segment_id]
                stable_point_idx = [0, len(self.output.reaction_string_images[segment_id]) - 1]
                saddle_point_idx: List[int] = []
    
                # To make mypy happy, declaration is needed.
                idx = self.output.transition_state_indices[segment_id]
                if idx is not None:
                    saddle_point_idx.append(idx)
            else:
                flatten_atoms, stable_point_idx, saddle_point_idx = self._flatten_images()
    
            energies = np.array([atoms.ase_atoms.get_potential_energy() for atoms in flatten_atoms])
            fig = go.Figure()
            steps = np.arange(len(energies))
            fig.add_trace(
                go.Scatter(
                    x=steps,
                    y=energies,
                    showlegend=False,
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=steps[stable_point_idx],
                    y=energies[stable_point_idx],
                    mode="markers",
                    name="stable point",
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=steps[saddle_point_idx],
                    y=energies[saddle_point_idx],
                    mode="markers",
                    name="saddle point",
                )
            )
    
            if plt_name is not None:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.to_gif)    def to_gif(self, gif_name: Optional[str] = None, rotation: str = "0x,0y,0z") -> None:
            """Creates gif file to show the movement of atoms along the reaction path.
    
            Args:
                gif_name (str or None, optional): The file name to which the reaction path will be
                    saved. Defaults to None.
                rotation (str, optional): The viewing angle. Defaults to '0x,0y,0z'.
            """
            fig = plt.figure(facecolor="white")
            ax = fig.add_subplot()
            ase.io.write(
                gif_name,
                [atoms.ase_atoms for path in self.output.reaction_string_images for atoms in path],
                format="gif",
                ax=ax,
                rotation=(rotation),
            )
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.to_traj)    def to_traj(self, traj_name: str) -> None:
            """Creates ase.io.Trajectory to save the reaction path.
    
            Args:
                traj_name (str): The file name to which the NEB trajectory will be saved.
            """
            traj = Trajectory(traj_name, "w")
            for path in self.output.reaction_string_images:
                for atoms in path:
                    traj.write(atoms.ase_atoms)
            traj.close()
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature)class ReactionStringFeature(FeatureBase):
        """Reliable string method that surely connect IS and FS.
    
        All local minima are true local minima. All local maxima is the result of dimer method.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/reaction_string.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature.__init__)    def __init__(
            self,
            optimize_is: bool = True,
            optimize_fs: bool = True,
            fmax_rp: float = 0.05,
            fmax_ts: float = 0.05,
            fmax_rd: float = 0.01,
            fmax_eq: float = 0.001,
            local_extrema_tol: float = 1e-3,
            interesting: Optional[Callable[[List[List[Atoms]]], List[bool]]] = None,
            timeout: Optional[timedelta] = None,  # TODO just float
            logfile_rp: Optional[str] = None,
            logfile_eq: Optional[str] = None,
            observers_rp: Optional[List[Callable[[], None]]] = None,
            observers_eq: Optional[List[Callable[[], None]]] = None,
            remove_rotation_and_translation: Optional[bool] = None,
            max_workers: Optional[int] = None,
            optimizer_eq: Union[OptimizerBuilder, str] = "FIRELBFGS",  # type: ignore
            optimizer_rp: Union[OptimizerBuilder, str] = "FIRES",  # type: ignore
            k: float = 0.1,
            dx_ts: float = 0.02,
            dx_rp: float = 0.3,
            dx_rm: float = 1.0,
            estimator_fn: Optional[EstimatorFnType] = None,
            dump_directory: Optional[Path] = None,
            dump_frequency: Optional[timedelta] = None,
        ) -> None:
            """Initialize an instance.
    
            Args:
                optimize_is (bool, optional):
                    Optimize initial structure. Defaults to True.
                optimize_fs (bool, optional):
                    Optimize final structure. Defaults to True.
                fmax_rp (float, optional):
                    Convergence criteria for reaction path. Defaults to 0.05.
                fmax_ts (float, optional):
                    Convergence criteria for transition states.
                    If math.isinf(fmax_ts), climbing algorithm will not be executed while optimizing.
                    Defaults to 0.05.
                fmax_rd (float, optional):
                    Convergence criterial for rate determining step (highest TS).
                    If math.isinf(fmax_rd), climbing algorithm will not be executed while optimizing.
                    Defaults to 0.01.
                fmax_eq (float, optional):
                    Convergence criteria for local minimas. Defaults to 0.001.
                local_extrema_tol (float, optional):
                    Local extrema smaller than this will be ignored. Defaults to 1e-3.
                interesting (list[list[Atoms]] -> list[bool] or None, optional):
                    Interesting path or not. Defaults to None.
                timeout (timedelta or None, optional):
                    Timeout. Defaults to None.
                logfile_rp (str or None, optional):
                    Logfile for reaction path finding algorithms.
                    If *logfile* is a string, a file with that name will be opened.
                    Use '-' for stdout. Defaults to None.
                logfile_eq (str or None, optional):
                    Logfile for local minima optimization.
                    If *logfile* is a string, a file with that name will be opened.
                    Use '-' for stdout. Defaults to None.
                observers_rp (list[ -> None] or None, optional):
                    Called in each step of reaction path optimzation. Defaults to None.
                observers_eq (list[ -> None] or None, optional):
                    Called after local minimas are optimized. Defaults to None.
                remove_rotation_and_translation (bool or None, optional):
                    True actives NEB-TR for removing translation and
                    rotation during NEB. By default applied non-periodic
                    systems Defaults to None.
                max_workers (int or None, optional):
                    Max number of concurrent threads. If it is not specified,
                    it is automatically determined based on the number of input atoms.
                optimizer_eq (OptimizerBuilder or str, optional):
                    Optimizer generator for structural optimization. Defaults to "FIRELBFGS".
                    Passing OptimizerBuilder type objects will be deprecated in a future version.
                optimizer_rp (OptimizerBuilder or str, optional):
                    Optimizer generator for reaction path optimization. Defaults to "FIRE".
                    Passing OptimizerBuilder type objects will be deprecated in a future version.
                k (float, optional):
                    NEB spring constant for String and NEB hybrid algorithm. Defaults to 0.1.
                dx_ts (float, optional):
                    Distance between images for climbing string. Defaults to 0.02.
                dx_rp (float, optional):
                    Distance between images for reaction path. Defaults to 0.3.
                dx_rm (float, optional):
                    Images closer than this value will be linearly interpolated to skip connection.
                    Defaults to 1.0.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail. Defaults to None.
                dump_directory (Path or None, optional):
                    Folder to save images during optimization.
                    If `dump_directory` is not `None`, you must specify `dump_frequency` as well.
                dump_frequency (timedelta or None, optional):
                    Interval to dump images on `dump_directory`.
                    If `dump_frequency` is not `None`, you must specify `dump_directory` as well.
            """
            super(ReactionStringFeature, self).__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            with self.init_scope():
                self.optimize_is = optimize_is
                self.optimize_fs = optimize_fs
                self.fmax_rp = fmax_rp
                self.fmax_ts = fmax_ts
                self.fmax_rd = fmax_rd
                self.fmax_eq = fmax_eq
                self.local_extrema_tol = local_extrema_tol
                self.interesting = interesting
                self.timeout = timeout
                self.logfile_rp = logfile_rp
                self.logfile_eq = logfile_eq
                self.observers_rp = observers_rp
                self.observers_eq = observers_eq
                self.remove_rotation_and_translation = remove_rotation_and_translation
                self.max_workers = max_workers
    
                # This refactoring is preferred for ease of save and reload.
                if isinstance(optimizer_eq, str):
                    self.optimizer_eq_str = optimizer_eq
                    if optimizer_eq == "FIRELBFGS":
                        self.optimizer_eq = FIRELBFGS
                    elif optimizer_eq == "FIRES":
                        self.optimizer_eq = FIRES
                    else:
                        raise ValueError("optimizer_eq must be 'FIRELBFGS' or 'FIRES'.")
                elif optimizer_eq in [FIRELBFGS, FIRES]:
                    warnings.warn(
                        "Use of ReactionString optimizer object for optimizer_eq is deprecated. "
                        'Please use a string, e.g. "FIRELBFGS", instead.'
                    )
                    self.optimizer_eq_str = optimizer_eq.__name__
                    self.optimizer_eq = optimizer_eq
                else:
                    raise TypeError(f"optimizer_eq must be str or Optimizer, got {type(optimizer_eq).__name__}")
    
                # This refactoring is preferred for ease of save and reload.
                if isinstance(optimizer_rp, str):
                    self.optimizer_rp_str = optimizer_rp
                    if optimizer_rp == "FIRELBFGS":
                        self.optimizer_rp = FIRELBFGS
                    elif optimizer_rp == "FIRES":
                        self.optimizer_rp = FIRES
                    else:
                        raise ValueError("optimizer_rp must be 'FIRELBFGS' or 'FIRES'.")
                elif optimizer_rp in [FIRELBFGS, FIRES]:
                    warnings.warn(
                        'Use of ReactionString optimizer object for optimizer_rp is deprecated. Please use a string, e.g. "FIRES", instead.'
                    )
                    self.optimizer_rp_str = optimizer_rp.__name__
                    self.optimizer_rp = optimizer_rp
                else:
                    raise TypeError(f"optimizer_rp must be str or Optimizer, got {type(optimizer_rp).__name__}")
    
                self.k = k
                self.dx_ts = dx_ts
                self.dx_rp = dx_rp
                self.dx_rm = dx_rm
                self.estimator_fn = estimator_fn
                self.dump_directory = dump_directory
                self.dump_frequency = dump_frequency
                self.dumpimages_enabled = dump_directory is not None or dump_frequency is not None
    
    
    
        def __call__(
            self,
            images: List[List[Union[MatlantisAtoms, ASEAtoms]]],
        ) -> ReactionStringFeatureResult:
            """Get the reaction path between the initial and end images.
    
            Args:
                images (list[list[MatlantisAtoms or ASEAtoms]]):
                    Images defining path from initial to final state.
            Returns:
                ReactionStringFeatureResult : The calculation result.
            """
            init_images = [[at.ase_atoms if isinstance(at, MatlantisAtoms) else at for at in atoms_list] for atoms_list in images]
    
            estimator_fn = self.estimator_fn
    
            def build_calculator(_: int) -> Calculator:
                return get_calculator(estimator_fn=estimator_fn)
    
            # NOTE: if ``max_workers`` is not specified,
            # it is automatically determined based on the number of input atoms.
            # s.t. max_workers in [1, 50]
            if self.max_workers is None:
                n_atoms = len(images[0][0])
                max_workers = max(1, min(50, 10000 // n_atoms))
            else:
                max_workers = self.max_workers
    
            if isinstance(self.observers_rp, list) and any(isinstance(observer, DumpImages) for observer in self.observers_rp):
                raise RuntimeError(
                    "You cannot provide `DumpImages` to `observers_rp` directly. Please specify"
                    " `dump_frequency` and `dump_directory when instantiating `ReactionStringFeature`."
                )
    
            if isinstance(self.observers_eq, list) and any(isinstance(observer, DumpImages) for observer in self.observers_eq):
                raise RuntimeError("`DumpImages` is not supported for `observers_eq`.")
    
            if self.dumpimages_enabled:
                if self.dump_directory is None:
                    raise RuntimeError("`dump_directory` must be provided if `dump_frequency` is specified.")
                if self.dump_frequency is None:
                    raise RuntimeError("`dump_frequency` must be provided if `dump_directory` is specified.")
                self.observers_rp = self.observers_rp or []
                self.observers_rp.append(
                    DumpImages(
                        images=init_images,
                        dump_directory=self.dump_directory,
                        dump_frequency=self.dump_frequency,
                    )
                )
    
            with ReactionPath(
                images=init_images,
                build_calculator=build_calculator,
                optimize_is=self.optimize_is,
                optimize_fs=self.optimize_fs,
                fmax_rp=self.fmax_rp,
                fmax_ts=self.fmax_ts,
                fmax_rd=self.fmax_rd,
                fmax_eq=self.fmax_eq,
                local_extrema_tol=self.local_extrema_tol,
                interesting=self.interesting,
                timeout=self.timeout,
                logfile_rp=self.logfile_rp,
                logfile_eq=self.logfile_eq,
                observers_rp=self.observers_rp,
                observers_eq=self.observers_eq,
                remove_rotation_and_translation=self.remove_rotation_and_translation,
                max_workers=max_workers,
                optimizer_eq=self.optimizer_eq,
                optimizer_rp=self.optimizer_rp,
                k=self.k,
                dx_ts=self.dx_ts,
                dx_rp=self.dx_rp,
                dx_rm=self.dx_rm,
            ) as rp:
                rp.run()
    
            if isinstance(self.observers_rp, list) and isinstance(self.observers_rp[-1], DumpImages):
                self.observers_rp.pop()
    
            reaction_path_images = rp.images
            reaction_string_images = self._convert_to_matlantis_atoms(reaction_path_images)
            output = ReactionStringFeatureOutput(
                num_segments=len(reaction_string_images),
                reaction_string_images=reaction_string_images,
                transition_state_indices=rp.ts,
                imaginary_eigenvectors=rp.imaginary_eigenvectors,
                converged=rp.converged,
            )
            result = ReactionStringFeatureResult(
                feature=deepcopy(self),
                call_params={"images": self._convert_to_matlantis_atoms(images)},
                output=output,
                calculator_info=get_calculator_info(estimator_fn=estimator_fn),
                version_info=get_version_info(),
                user_info=get_user_info(),
            )
            return result
    
        def _convert_to_matlantis_atoms(
            self,
            reaction_path_images: List[List[Union[Atoms, MatlantisAtoms]]],
        ) -> List[List[MatlantisAtoms]]:
            """Convert ASE Atoms to MatlantisAtoms.
    
            Atoms in `reaction_path_images` might not have a calculator,
            typically when `timeout` is extremely short and reactionstring
            finishes their calculation before attatching that. To avoid
            raising an error in such case, we manually attach a calculator
            if Atoms object doesn't have a calculator.
    
            Args:
                reaction_path_images (list[list[Atoms or MatlantisAtoms]]):
                    Images defining path from initial to final state.
            Returns:
                list[list[MatlantisAtoms]] : The calculation result of `MatlantisAtoms`.
            """
            ret = []
            for path in reaction_path_images:
                path_matlantis_atoms = []
                for atoms in path:
                    if isinstance(atoms, Atoms):
                        if atoms.calc is None:
                            warnings.warn("calculator is not attached.")
                            atoms.calc = get_calculator(self.estimator_fn)
                        path_matlantis_atoms.append(MatlantisAtoms(atoms))
                    else:
                        atoms.ase_atoms.calc = get_calculator(self.estimator_fn)
                        path_matlantis_atoms.append(atoms)
                ret.append(path_matlantis_atoms)
    
            return ret
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the ReactionStringFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized ReactionStringFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "optimize_is": self.optimize_is,
                "optimize_fs": self.optimize_fs,
                "fmax_rp": self.fmax_rp,
                "fmax_ts": self.fmax_ts,
                "fmax_rd": self.fmax_rd,
                "fmax_eq": self.fmax_eq,
                "local_extrema_tol": self.local_extrema_tol,
                "interesting": MatlantisNonSerializable.from_param("interesting", self.interesting),
                "timeout": MatlantisNonSerializable.from_param("timeout", self.timeout),
                "logfile_rp": self.logfile_rp,
                "logfile_eq": self.logfile_eq,
                "observers_rp": MatlantisNonSerializable.from_param("observers_rp", self.observers_rp),
                "observers_eq": MatlantisNonSerializable.from_param("observers_eq", self.observers_eq),
                "remove_rotation_and_translation": self.remove_rotation_and_translation,
                "max_workers": self.max_workers,
                "optimizer_eq": self.optimizer_eq_str,
                "optimizer_rp": self.optimizer_rp_str,
                "k": self.k,
                "dx_ts": self.dx_ts,
                "dx_rp": self.dx_rp,
                "dx_rm": self.dx_rm,
                "estimator_fn": MatlantisNonSerializable.from_param("estimator_fn", self.estimator_fn),
                "dump_directory": self.dump_directory,
                "dump_frequency": self.dump_frequency,
            }
            return to_serialized_dict(res)
    
    
    
