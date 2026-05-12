# Source code for matlantis_features.features.phonon.mode
    
    
    import os
    import warnings
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional, Tuple, Union
    
    import ase
    import matplotlib.pyplot as plt
    import numpy as np
    from ase import Atoms
    from ase.io import Trajectory
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.force_constant import ForceConstantFeatureResult
    from matlantis_features.features.phonon.utils import PhononFrequency
    from matlantis_features.utils import FeatureCost
    
    # from ase.io.animation import write_gif
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html#matlantis_features.features.phonon.mode.PostPhononModeFeatureResult)@dataclass
    class PostPhononModeFeatureResult:
        """A dataclass for result of PostPhononModeFeature."""
    
        k_point: Tuple[float, float, float]
        repeat_of_cell: Tuple[int, int, int]
        amplitude: float
        n_images: int
        mode_indices: List[int]
        frequencies: List[float]
        trajectories: List[List[Atoms]]
        q_point: Optional[Tuple[float, float, float]] = None
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html#matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.plot)    def plot(self, save_path: str, rotation: str = "0x,0y,0z") -> None:
            """Creates gif files to show the movement of atoms in each phonon mode.
    
            Args:
                save_path (str): The directory to which the gif files will be saved. The file will be
                    named as 'mode' + band index + '.gif'. The band index is given in ascending order of
                    phonon frequency.
                rotation (str, optional): The viewing angle. Defaults to '0x,0y,0z'.
            """
            os.makedirs(save_path, exist_ok=True)
            for i, t in zip(self.mode_indices, self.trajectories):
                fig = plt.figure(facecolor="white")
                ax = fig.add_subplot()
                ase.io.write(f"{save_path}/mode_{i}.gif", t, format="gif", ax=ax, rotation=(rotation))
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html#matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.to_traj)    def to_traj(self, save_path: str) -> None:
            """Create ase.io.Trajectory to save the movement of atoms in each phonon mode.
    
            Args:
                save_path (str): The directory to which the gif files will be saved. The file will be
                    named as 'mode' + band index + '.traj'. The band index is given in ascending order of
                    phonon frequency.
            """
            os.makedirs(save_path, exist_ok=True)
            for i, t in zip(self.mode_indices, self.trajectories):
                traj = Trajectory(f"{save_path}/mode_{i}.traj", "w")
                for atoms in t:
                    traj.write(atoms)
                traj.close()
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the PostPhononModeFeatureResult.
    
            Returns:
                dict[str, Any] : Phonon modes in dictionary format.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.mode.PostPhononModeFeature.html#matlantis_features.features.phonon.mode.PostPhononModeFeature)class PostPhononModeFeature(FeatureBase):
        """The matlantis-feature to visualize how atoms move under a certain phonon mode."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.mode.PostPhononModeFeature.html#matlantis_features.features.phonon.mode.PostPhononModeFeature.__init__)    def __init__(self) -> None:
            """Initialize an instance."""
            super(PostPhononModeFeature, self).__init__()
    
    
    
        def __call__(
            self,
            force_constant: ForceConstantFeatureResult,
            k_point: Tuple[float, float, float] = (0.0, 0.0, 0.0),
            repeat_of_cell: Tuple[int, int, int] = (3, 3, 3),
            mode_indices: Optional[Union[int, List[int]]] = None,
            n_images: int = 20,
            amplitude: float = 5.0,
            q_point: Optional[Tuple[float, float, float]] = None,
        ) -> PostPhononModeFeatureResult:
            """Calculates the phonon mode from the force constant.
    
            Args:
                force_constant (ForceConstantFeatureResult): The force contant result obtained from
                    matlantis_features.features.phonon.RunForceConstantFeature
                k_point (tuple[float, float, float], optional): The k-point in which the phonon mode
                    will be calculated and the atom movement will be viewed. Defaults to [0., 0., 0.].
                repeat_of_cell (tuple[int, int, int], optional): To view the phonon mode at a large
                    supercell. Defaults to [3, 3, 3].
                mode_indices (int or list[int] or None, optional): The index of branches.
                    If None is provided, the phonon modes will be calculated in all branches.
                n_images (int, optional): Number of images in one period. Defaults to 20.
                amplitude (float, optional):
                    The amplitude of vibration. Defaults to 1.0.
                q_point (tuple[float, float, float] or None, optional): Deprecated parameter.
                    The q-vector in which the phonon mode will be calculated and the atom movement
                    will be viewed. Defaults to [0., 0., 0.].
    
            Returns:
                PostPhononModeFeatureResult:
                    A set of structure images that show the periodic movement of each atom under phonon
                    modes of given a k-point.
            """
            if q_point is not None:
                k_point = q_point
                warnings.warn(
                    "The parameter 'q_point' is deprecated. Please use 'k_point' instead.",
                    FutureWarning,
                )
    
            assert len(k_point) == 3
            assert len(repeat_of_cell) == 3
    
            pf = PhononFrequency(
                force_constant=force_constant.force_constant,
                supercell=force_constant.supercell,
                unit_cell_atoms=force_constant.unit_cell_atoms,
            )
            n_atoms = pf.n_atoms
            # TODO check the q in brillouin zone
            frequency, mode = pf._frequency_mode_q(np.array(k_point), return_mode=True)  # type: ignore
    
            atoms = force_constant.unit_cell_atoms.ase_atoms.copy() * repeat_of_cell
    
            n_cell = np.prod(repeat_of_cell)
            positions = atoms.get_positions()
            R = np.indices(repeat_of_cell).reshape(3, -1)
    
            if mode_indices is None:
                mode_indices = list(range(n_atoms * 3))
            elif isinstance(mode_indices, int):
                mode_indices = [mode_indices]
    
            traj: List[List[Atoms]] = [[] for _ in mode_indices]
    
            for i in mode_indices:
                displacement = np.tile(mode[i].real, [n_cell, 1]) * amplitude
    
                for x in np.linspace(0, 2 * np.pi, n_images, endpoint=False):
                    phase = np.exp(-2.0j * np.pi * np.dot(k_point, R) + 1.0j * x).real.repeat(n_atoms)
                    positions_disp = positions + displacement * phase[:, np.newaxis]
    
                    atoms_disp = atoms.copy()
                    atoms_disp.set_positions(positions_disp)
                    traj[i].append(atoms_disp)
    
            results = PostPhononModeFeatureResult(
                k_point=k_point,
                repeat_of_cell=repeat_of_cell,
                amplitude=amplitude,
                n_images=n_images,
                mode_indices=mode_indices,
                frequencies=frequency[mode_indices].tolist(),
                trajectories=traj,
                q_point=k_point,
            )
    
            return results
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): The number of atoms.
            Returns:
                FeatureCost :
                    A dataclass containing cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
