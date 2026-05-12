# Source code for matlantis_features.features.md.post_features.diffusion
    
    
    import warnings
    from collections import OrderedDict
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional, Tuple, Union
    
    import numpy as np
    import plotly
    import plotly.graph_objects as go
    from ase.data import chemical_symbols
    from ase.md.analysis import DiffusionCoefficient
    from ase.units import fs
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md.md import MDFeatureResult
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import (
        diffusion_coefficient_converter,
        diffusion_coefficient_converter_old,
    )
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult)@dataclass
    class PostMDDiffusionFeatureResult:
        """A dataclass for result of PostMDDiffusionFeature."""
    
        diffusion_coefficient: Dict[str, Dict[str, float]]
        diffusion_coefficient_std: Dict[str, Dict[str, float]]
        diffusion_coefficient_molecule: Dict[str, Dict[str, float]]
        diffusion_coefficient_molecule_std: Dict[str, Dict[str, float]]
        mean_squared_displacement: Dict[str, np.ndarray]
        mean_squared_displacement_molecule: Dict[str, np.ndarray]
        timestep: float
        init_time: float
        stride: int
        atom_indices: Optional[List[int]]
        molecule: Union[bool, List[List[int]]]
        number_of_segments: int
        direction: Optional[np.ndarray]
        method: str
        effective_msd_range: Optional[Tuple[float, float]]
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.plot)    def plot(
            self,
            plt_name: Optional[str] = None,
            molecule: bool = False,
            show_fit_line: bool = False,
            show_effective_range: bool = False,
        ) -> go.Figure:
            """Plot the mean square displacement.
    
            Args:
                plt_name (str or None, optional):
                  File name to write the plot. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                  'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                  Defaults to None.
                molecule (bool, optional):
                  Plot the mean sqare displacement of molecular center of mass. Defaults to None.
                show_fit_line (bool, optional):
                  Plot the fitting line. Defaults to False.
                show_effective_range (bool, optional):
                  Plot the effective range of msd that used for the linear fitting. Defaults to False.
            Returns:
                go.Figure : Resulting plotly's graph object.
            """
            if molecule:
                if len(self.mean_squared_displacement_molecule) == 0:
                    raise ValueError("No molecular information is found in the calculation result.")
                msd = self.mean_squared_displacement_molecule
            else:
                msd = self.mean_squared_displacement
    
            n_steps = msd.values().__iter__().__next__().shape[0]
            time = np.arange(n_steps) * self.timestep
    
            init, end = _get_init_end_step(n_steps, effective_range=self.effective_msd_range)
    
            fig = go.Figure()
            colors = plotly.colors.qualitative.Plotly
            for i, (k, m) in enumerate(msd.items()):
                y = m[:, :3].sum(axis=1)
                fig.add_trace(go.Scatter(x=time, y=y, name=k, line=dict(color=colors[i])))
    
                if show_fit_line:
                    f = np.poly1d(np.polyfit(time[init:end], y[init:end], deg=1))
                    fig.add_trace(
                        go.Scatter(
                            x=time[[init, end]],
                            y=f(time[[init, end]]),
                            line=dict(width=4, color=colors[i], dash="dash"),
                            showlegend=False,
                        )
                    )
    
                if show_effective_range:
                    fig.add_vrect(
                        x0=time[init],
                        x1=time[end],
                        fillcolor="gray",
                        opacity=0.1,
                        annotation_text="effective range",
                        annotation_position="top left",
                    )
    
            fig.update_xaxes(title="time [fs]", rangemode="tozero")
            fig.update_yaxes(title=r"$\textrm{Mean squared displacement} [A^2] $", rangemode="tozero")
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Convert the result to the dict object.
    
            Returns:
                dict[str, Any] : Resulting dict object.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature)class PostMDDiffusionFeature(FeatureBase):
        """The matlantis-feature for calculating the diffusion coefficient using the MD trajectory.
    
        .. note::
            The breaking change is made in ``PostMDDiffusionFeature`` at version 0.8.0.
            The new implementation is fast,
            and enables a more comprehensive analysis of diffusion properties.
            The old version is renamed as
            :class:`~matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature`.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/diffusion.html>`_.
        """
    
        def __call__(
            self,
            md_results: MDFeatureResult,
            init_time: float = 0.0,
            stride: int = 1,
            atom_indices: Optional[List[int]] = None,
            molecule: Union[bool, List[List[int]]] = False,
            number_of_segments: int = 1,
            direction: Optional[np.ndarray] = None,
            method: str = "normal",
            effective_msd_range: Optional[Tuple[float, float]] = None,
        ) -> PostMDDiffusionFeatureResult:
            """Calculate the diffusion coefficient of each element from the MD trajectory.
    
            Args:
                md_results (MDFeatureResult):
                  MD feature's result object containing the trajectory.
                init_time (float, optional):
                  Trajectory frames before this time (in fs) is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                stride (int, optional):
                  Size of the skip strie of the trajectory frames.
                  Defaults to 1 (i.e., all frames are used for calculation).
                atom_indices (list[int] or None, optional):
                  The indices of atoms whose diffusion coefficient is to be calculated explicitly.
                  Defaults to None.
                molecule (bool or list[list[int]], optional):
                  The indices of atoms who form a molecule. If it is used, the diffusion coefficient
                  of the center of mass of molecules will be calculated. The results will be saved in
                  the `PostMDDiffusionFeatureResult.diffusion_coefficient_molecule`. If True, all the
                  atoms will be regard as one molecule. Defaults to False.
                number_of_segments (int, optional):
                  Divides the given trajectory in to segments to allow statistical analysis, such as
                  standard deviation of diffuion coefficient.
                  Defaults to 1.
                direction (np.ndarray or None, optional):
                  Calculate the diffusion coefficient along specific directions. It must be a nx3 array
                  where n is number of directions. The result will be saved in
                  `PostMDDiffusionFeatureResult.diffusion_coefficient` as, e.g. `H_direc_1`,
                  `O_direc_2` etc. If None, only the diffusion coefficient along x, y and z directions
                  is calculated. Defaults to None
                method (str, optional):
                  The method to calculate mean square displacement (MSD). Now, two methods, i.e.
                  "normal" and "segment", are supported. Defaults to "normal".
                effective_msd_range (tuple[float, float] or None, optional):
                  A tuple of two floats representing the range of MSD sequence
                  to be used in the calculation of diffusion coefficient.
                  The values should be in the range of [0.0, 1.0].
                  If None, all the mean square displacements will be used.
                  Defaults to None.
    
            Returns:
                PostMDDiffusionFeatureResult : Result dataclass object.
            """
            assert init_time >= 0.0
            assert stride >= 1
            assert method in ["normal", "segment"]
            if effective_msd_range is not None:
                assert len(effective_msd_range) == 2
                assert effective_msd_range[0] >= 0.0 and effective_msd_range[0] <= 1.0
                assert effective_msd_range[1] >= 0.0 and effective_msd_range[1] <= 1.0
                assert effective_msd_range[0] < effective_msd_range[1]
    
            traj_raw = md_results.get_traj_obj(self.get_savedir_from_ctx())
            try:
                time_raw = np.array([atoms.info["time"] for atoms in traj_raw])
            except KeyError:
                raise ValueError("The trajectory does not contains time information.")
    
            if not np.allclose(np.diff(time_raw), np.mean(np.diff(time_raw))):
                raise ValueError("The time interval is not all the same in this trajectory.")
    
            if time_raw.max() <= init_time:
                raise ValueError("init_time is larger than the length of trajectory.")
    
            start_step = np.where(time_raw >= init_time)[0][0]
    
            traj = traj_raw[start_step::stride]
            time = time_raw[start_step::stride]
            time -= time[0]
            time_interval = np.mean(np.diff(time))
    
            n_steps = len(traj)
            n_atoms = len(traj[0])
            positions_array = np.stack([atoms.get_positions() for atoms in traj])
            numbers = traj[0].get_atomic_numbers()
            masses = traj[0].get_masses().reshape([1, -1, 1])
            n_steps_per_segment = n_steps // number_of_segments
    
            if atom_indices:
                if max(atom_indices) >= n_atoms:
                    raise ValueError(f"The index {max(atom_indices)} is larger than the total number of atoms {n_atoms} in 'atom_indices'.")
                elif min(atom_indices) < 0:
                    raise ValueError("The index should be positive integer in 'atom_indices'.")
    
            molecule_indices = self._process_molecule(molecule, n_atoms)
    
            element_indices = _get_element_indices(numbers, atom_indices)
            msd_seg = []
            msd_mol_seg = []
            for i in range(number_of_segments):
                msd_, msd_mol_ = get_msd(
                    positions_array=positions_array[n_steps_per_segment * i : n_steps_per_segment * (i + 1)],
                    masses=masses,
                    element_indices=element_indices,
                    com_indices=molecule_indices,
                    direction=direction,
                    method=method,
                )
                msd_seg.append(msd_)
                msd_mol_seg.append(msd_mol_)
    
            msd = np.mean(msd_seg, axis=0)
            msd_mol = np.mean(msd_mol_seg, axis=0)
    
            mean_squared_displacement = {k: v for k, v in zip(element_indices.keys(), msd.transpose([1, 0, 2]))}
            mean_squared_displacement_molecule = {f"molecule_{i}": v for i, v in enumerate(msd_mol.transpose([1, 0, 2]))}
    
            diffusion_coefficient, diffusion_coefficient_std = _get_diffusion_coefficients(
                n_steps_per_segment, time_interval, mean_squared_displacement, effective_msd_range
            )
    
            (
                diffusion_coefficient_molecule,
                diffusion_coefficient_molecule_std,
            ) = _get_diffusion_coefficients(
                n_steps_per_segment,
                time_interval,
                mean_squared_displacement_molecule,
                effective_msd_range,
            )
    
            for atom_type in diffusion_coefficient["A^2/fs"].keys():
                if diffusion_coefficient["A^2/fs"][atom_type] < 0.0:
                    warnings.warn(
                        "Negative diffusion coefficient detected",
                    )
                if diffusion_coefficient_std["A^2/fs"][atom_type] / diffusion_coefficient["A^2/fs"][atom_type] > 0.5:
                    warnings.warn(
                        "Large variance detected in fitting.",
                    )
    
            results = PostMDDiffusionFeatureResult(
                diffusion_coefficient=diffusion_coefficient,
                diffusion_coefficient_std=diffusion_coefficient_std,
                diffusion_coefficient_molecule=diffusion_coefficient_molecule,
                diffusion_coefficient_molecule_std=diffusion_coefficient_molecule_std,
                mean_squared_displacement=mean_squared_displacement,
                mean_squared_displacement_molecule=mean_squared_displacement_molecule,
                timestep=time_interval,
                init_time=init_time,
                stride=stride,
                atom_indices=atom_indices,
                molecule=molecule,
                number_of_segments=number_of_segments,
                direction=direction,
                method=method,
                effective_msd_range=effective_msd_range,
            )
    
            return results
    
        def _process_molecule(self, molecule: Union[bool, List[List[int]]], n_atoms: int) -> Optional[List[List[int]]]:
            molecule_indices: Optional[List[List[int]]]
            if isinstance(molecule, bool):
                if molecule is True:
                    warnings.warn(
                        "If 'molecule=True', all the atoms will be regarded as one molecule. "
                        "If it is not true, please specify the indices of atoms that form one molecule."
                        "For example, 'molecule=[[0,1,2], [3,4,5]]' indicates that the atom 0, 1 and 2 "
                        "form one molecule, and the atom 2, 3 and 4 form another molecule."
                    )
                    molecule_indices = [[i for i in range(n_atoms)]]
                elif molecule is False:
                    molecule_indices = None
            else:
                molecule_indices = molecule
                max_ind = max(sum(molecule_indices, []))
                if max_ind >= n_atoms:
                    raise ValueError(f"The index {max_ind} is larger than the total number of atoms {n_atoms} in 'molecule_indices'.")
                min_ind = min(sum(molecule_indices, []))
                if min_ind < 0:
                    raise ValueError("The index should be positive integer in 'molecule_indices'.")
            return molecule_indices
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
    
    def _get_diffusion_coefficients(
        n_steps: int,
        time_interval: float,
        mean_squared_displacement: Dict[str, np.ndarray],
        effective_msd_range: Optional[Tuple[float, float]] = None,
    ) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
        time = np.arange(n_steps) * time_interval
        init, end = _get_init_end_step(n_steps, effective_msd_range)
    
        units = ["A^2/fs", "cm^2/s"]
        diffusion_coefficient: Dict[str, Dict[str, float]] = {k: {} for k in units}
        diffusion_coefficient_std: Dict[str, Dict[str, float]] = {k: {} for k in units}
    
        for element, msd in mean_squared_displacement.items():
            m = msd[:, :3].sum(axis=1)
            D, D_std = fit_diffusion_coefficients(time[init:end], m[init:end], dimensions=3)
            D_unit = diffusion_coefficient_converter(D)
            D_std_unit = diffusion_coefficient_converter(D_std)
            for u in units:
                diffusion_coefficient[u][element] = D_unit[u]  # type: ignore
                diffusion_coefficient_std[u][element] = D_std_unit[u]  # type: ignore
            for i, m in enumerate(msd.T):
                if i < 3:
                    key = f"{element}_{'xyz'[i]}"
                else:
                    key = f"{element}_direc_{i - 3}"
                D, D_std = fit_diffusion_coefficients(time[init:end], m[init:end])
                D_unit = diffusion_coefficient_converter(D)
                D_std_unit = diffusion_coefficient_converter(D_std)
                for u in units:
                    diffusion_coefficient[u][key] = D_unit[u]  # type: ignore
                    diffusion_coefficient_std[u][key] = D_std_unit[u]  # type: ignore
        return diffusion_coefficient, diffusion_coefficient_std
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients.html#matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients)def fit_diffusion_coefficients(time: np.ndarray, msd: np.ndarray, dimensions: int = 1) -> Tuple[float, float]:
        """Obtain the diffusion coefficient from fitting the time and the mean squared displacement.
    
        Args:
            time (np.ndarray): The time of each step in fs.
            msd (np.ndarray): The mean squared displacement.
            dimensions (int, optional): The dimensions. Defaults to 1.
    
        Returns:
            Tuple[float, float]: The diffusion coefficient and its standard deviation.
        """
        p, cov = np.polyfit(time, msd, deg=1, cov=True)
        D = p[0] / (2 * dimensions)
        D_std = np.sqrt(cov[0, 0]) / (2 * dimensions)
        return D, D_std
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.get_msd.html#matlantis_features.features.md.post_features.diffusion.get_msd)def get_msd(
        positions_array: np.ndarray,
        masses: np.ndarray,
        element_indices: OrderedDict,
        com_indices: Optional[List[List[int]]] = None,
        direction: Optional[np.ndarray] = None,
        method: str = "normal",
    ) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        """Calculate the mean squared displacement.
    
        Args:
            positions_array (np.ndarray): The atomic position of each atom in each step.
            masses (np.ndarray): The atomic number of each atom.
            element_indices (OrderedDict): The indices of atoms of each element.
            com_indices (Optional[List[List[int]]], optional):
              The indices of atoms who form a molecule. If it is used, the diffusion coefficient
              of the center of mass of molecules will be calculated. Defaults to None.
            direction (Optional[np.ndarray], optional):
              Calculate the diffusion coefficient along specific directions. It must be a nx3 array
              where n is number of directions. If None, only the diffusion coefficient along
              x, y and z directions is calculated. Defaults to None
            method (str, optional):
              The method to calculate mean square displacement (MSD). Now, two methods, i.e.
              "normal" and "segment", are supported. Defaults to "normal".
        Returns:
            tuple[np.ndarray, np.ndarray or None]: the mean square displacement of each atom, the mean square
            displacement of each molecule.
        """
        steps, n_atoms = positions_array.shape[:2]
    
        if direction is not None:
            assert direction.shape[1] == 3
            direction_all = np.vstack((np.eye(3), direction))
            direction_norm = direction_all / np.linalg.norm(direction_all, axis=1, keepdims=True)
            positions_proj = positions_array.dot(direction_norm.T)
        else:
            positions_proj = positions_array
        n_direction = positions_proj.shape[2]
    
        msd_atoms = np.zeros([steps, n_atoms, n_direction])
        for i in range(1, steps):
            if method == "segment":
                disp_2 = (positions_proj[i:] - positions_proj[:-i]) ** 2
                msd_atoms[i] = np.mean(disp_2, axis=0)
            elif method == "normal":
                disp_2 = (positions_proj[i] - positions_proj[0]) ** 2
                msd_atoms[i] = disp_2
            else:
                raise ValueError("The 'segment' and 'normal' methods are supported!")
    
        msd_element = np.zeros([steps, len(element_indices), n_direction])
        for i, indices in enumerate(element_indices.values()):
            msd_element[:, i] = msd_atoms[:, indices].mean(axis=1)
    
        msd_com: Optional[np.ndarray]
        if com_indices:
            msd_com = np.zeros([steps, len(com_indices), n_direction])
            com_array = np.stack(
                [np.sum(positions_array[:, group] * masses[:, group], axis=1) / np.sum(masses[:, group]) for group in com_indices],
                axis=1,
            )
            if direction is not None:
                com_proj = com_array.dot(direction_norm.T)
            else:
                com_proj = com_array
            for i in range(1, steps):
                if method == "segment":
                    disp_2 = (com_proj[i:] - com_proj[:-i]) ** 2
                    msd_com[i] = np.mean(disp_2, axis=0)  # type: ignore
                elif method == "normal":
                    disp_2 = (com_proj[i] - com_proj[0]) ** 2
                    msd_com[i] = disp_2  # type: ignore
                else:
                    raise ValueError("The 'segment' and 'normal' methods are supported!")
        else:
            msd_com = np.empty([steps, 0, n_direction])
        return msd_element, msd_com
    
    
    
    
    def _get_element_indices(numbers: np.ndarray, atom_indices: Optional[List[int]] = None) -> OrderedDict:
        element_indices = OrderedDict()
        for n in np.unique(numbers):
            ind = np.where(numbers == n)[0].tolist()
            if atom_indices:
                ind = _intersection(ind, atom_indices)
            if len(ind) > 0:
                element_indices[chemical_symbols[n]] = ind
        return element_indices
    
    
    def _get_init_end_step(n_steps: int, effective_range: Optional[Tuple[float, float]] = None) -> Tuple[int, int]:
        if effective_range:
            assert effective_range[0] >= 0
            assert effective_range[1] <= 1.0
            assert effective_range[1] > effective_range[0]
            return int(effective_range[0] * n_steps), int(effective_range[1] * n_steps) - 1
        else:
            return 0, n_steps - 1
    
    
    def _intersection(list1: List[int], list2: List[int]) -> List[int]:
        return list(set(list1).intersection(list2))
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult)@dataclass
    class OldPostMDDiffusionFeatureResult:
        """A dataclass for result of PostMDDiffusionFeature."""
    
        diffusion_coefficient: Dict[str, Dict[str, float]]
        diffusion_coefficient_std: Dict[str, Dict[str, float]]
        mean_squared_displacement: Dict[str, np.ndarray]
        timestep: float
        init_time: float
        stride: int
        number_of_segments: int
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plot the result.
    
            Args:
                plt_name (str or None, optional):
                  File name to write the plot. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                  'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                  Defaults to None.
            Returns:
                go.Figure : Resulting plotly's graph object.
            """
            msd = self.mean_squared_displacement
            n_step = msd.values().__iter__().__next__().shape[1]
            time = np.arange(n_step) * self.timestep
            fig = go.Figure()
            colors = plotly.colors.qualitative.Plotly
            for i, (k, v) in enumerate(msd.items()):
                mean = v.mean(0)
                std = v.std(0)
                fig.add_trace(go.Scatter(x=time, y=mean, name=k, marker=dict(color=colors[i])))
                go.Scatter(
                    x=time + time[::-1],
                    y=(mean + 5 * std) + (mean - 5 * std)[::-1],
                    fill="toself",
                    fillcolor=colors[i],
                    line=dict(color=colors[i]),
                    hoverinfo="skip",
                    showlegend=False,
                )
            fig.update_xaxes(title="time [fs]")
            fig.update_yaxes(title="Mean squared displacement [A^2]")
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Convert the result to the dict object.
    
            Returns:
                dict[str, Any] : Resulting dict object.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature)class OldPostMDDiffusionFeature(FeatureBase):
        """The matlantis-feature for calculating the diffusion coefficient using the MD trajectory."""
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.__init__)    def __init__(self) -> None:
            super(OldPostMDDiffusionFeature, self).__init__()
            warnings.warn(
                "The OldPostMDDiffusionFeature will be deprecated in the future. Please use PostMDDiffusionFeature instead.",
                FutureWarning,
            )
    
    
    
        def __call__(
            self,
            md_results: MDFeatureResult,
            init_time: float = 0.0,
            stride: int = 1,
            atom_indices: Optional[List[int]] = None,
            molecule: bool = False,
            number_of_segments: int = 1,
        ) -> OldPostMDDiffusionFeatureResult:
            """Calculate the diffusion coefficient from the MD trajectory. \
            This feature class use the ASE's DiffusionCoefficent implementation.
    
            Args:
                md_results (MDFeatureResult):
                  MD feature's result object containing the trajectory.
                init_time (float, optional):
                  Trajectory frames before this time (in fs) is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                stride (int, optional):
                  Size of the skip strie of the trajectory frames.
                  Defaults to 1 (i.e., all frames are used for calculation).
                atom_indices (list[int] or None, optional):
                  The indices of atoms whose diffusion coefficient is
                  to be calculated explicitly.
                  Defaults to None.
                molecule (bool, optional):
                  Indicate if we are studying a molecule instead of atoms,
                  therefore use centre of mass in calculations.
                  Defaults to False.
                number_of_segments (int, optional):
                  Divides the given trajectory in to segments to allow statistical analysis.
                  Defaults to 1.
            Returns:
                OldPostMDDiffusionFeatureResult : Result dataclass object.
            """
            assert init_time >= 0.0
            assert stride >= 1
            assert number_of_segments >= 1
    
            traj = md_results.get_traj_obj(self.get_savedir_from_ctx())
            n_atoms = len(traj[0])
            if atom_indices:
                if max(atom_indices) >= n_atoms:
                    raise IndexError(f"The index {max(atom_indices)} is larger than the total number of atoms {n_atoms}.")
                elif min(atom_indices) < 0:
                    raise IndexError("The index should be positive integer.")
    
            timestep = traj[1].info["time"] - traj[0].info["time"]
            if traj[-1].info["time"] < init_time:
                raise ValueError("init_time is larger than the length of trajectory.")
    
            for i, atoms in enumerate(traj):
                if atoms.info["time"] >= init_time:
                    st = i
                    break
    
            sliced_traj = traj[st::stride]
    
            if len(sliced_traj) // number_of_segments <= 1:
                raise ValueError(
                    f"Length of sliced trajectory is {len(sliced_traj)}. "
                    "Please check the setting of 'init_time', 'stride' and 'number_of_segments'."
                )
            elif len(sliced_traj) // number_of_segments < 10:
                warnings.warn(f"Sliced trajectory is small ({len(sliced_traj)} snapshots). Please check if this is intended.")
    
            if molecule and (atom_indices is None):
                warnings.warn(
                    "Usually the atom_indices is needed to calculated the diffusion coefficient of molecules."
                    "Otherwise, it will calculate the displacement of the mass of center of the whole system."
                    "Please provide the list of atom indice that belong to one molecule in 'atom_indices'."
                )
    
            ase_dc = DiffusionCoefficient(
                traj=sliced_traj,
                timestep=timestep * stride * fs,
                atom_indices=atom_indices,
                molecule=molecule,
            )
            ase_dc.calculate(ignore_n_images=0, number_of_segments=number_of_segments)
            _slopes, _std = ase_dc.get_diffusion_coefficients()
            xyz = ase_dc.xyz_segment_ensemble_average
    
            for i, j in zip(_slopes, _std):
                if i < 0:
                    warnings.warn(
                        "Negative diffusion coefficient detected",
                    )
                if j / i > 0.5:
                    warnings.warn(
                        "Large variance detected in fitting.",
                    )
    
            slopes = diffusion_coefficient_converter_old(_slopes)
            std = diffusion_coefficient_converter_old(_std)
    
            slopes_dict, std_dict = {}, {}
            for k in slopes.keys():
                slopes_dict[k] = {atom_type: value for atom_type, value in zip(ase_dc.types_of_atoms, slopes[k])}
            for k in std.keys():
                std_dict[k] = {atom_type: value for atom_type, value in zip(ase_dc.types_of_atoms, std[k])}
    
            msd = {}
            for i, k in enumerate(ase_dc.types_of_atoms):
                msd[k] = xyz[:, i].sum(1) * 2.0
    
            results = OldPostMDDiffusionFeatureResult(
                diffusion_coefficient=slopes_dict,
                diffusion_coefficient_std=std_dict,
                mean_squared_displacement=msd,
                timestep=timestep * stride,
                init_time=init_time,
                stride=stride,
                number_of_segments=number_of_segments,
            )
    
            return results
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
