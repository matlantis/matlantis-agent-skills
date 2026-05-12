# Source code for matlantis_features.features.common.vibration
    
    
    from dataclasses import asdict, dataclass
    from math import sqrt
    from typing import Any, Dict, List, Optional, Tuple, Union
    
    import ase.units as units
    import numpy as np
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.common.hessian import _Hessian
    from matlantis_features.utils import FeatureCost, get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.VibrationFeatureResult.html#matlantis_features.features.common.vibration.VibrationFeatureResult)@dataclass
    class VibrationFeatureResult:
        """A dataclass for result of VibrationFeature."""
    
        atoms: MatlantisAtoms
        indices: Optional[List[int]]
        delta: float
        nfree: int
        method: str
        direction: str
        all_forces: np.ndarray
        force_constant: np.ndarray
        frequencies: np.ndarray
        image_frequencies: np.ndarray
        mode: np.ndarray
        zero_point_energy: float
        potential_energy: float
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the VibrationFeatureResult.
    
            Returns:
                dict[str, Any] : The vibration calculation results in dictionary format.
            """
            return asdict(self)
    
    
    
    
    _VIBRATION_DOC_STRING = """Initialize an instance.
    
            Args:
                delta (float, optional): The distance of finite displacement applied to each atom.
                    The unit is Angstrom. Defaults to 0.01.
                nfree (int, optional): The number of displacements per atom. Must be 2 or 4. Defaults
                    to 2.
                method (str, optional): The method to calculate the force constant from the displaced
                    supercells and their forces. This parameter can be 'Standard' or 'Frederiksen'.
                    Defaults to 'standard'.
                direction (str, optional): This parameter discribe how the displacment is found. Must
                    be  When 'direction=forward'/'direction=backward', each atom is displaced along the
                    positive/negative direction of the axis. When 'direction=center', each atom is
                    displaced along both positive and negative directions of the axis. Defaults to
                    'central'.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
                num_threads (int, optional): Number of threads for parallel calculation. Defaults to 10.
            """
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.VibrationFeature.html#matlantis_features.features.common.vibration.VibrationFeature)class VibrationFeature(FeatureBase):
        """The matlantis-feature for calculating molecular vibration.
    
        The vibration properties, including the force constant, vibration frequencies and \
        vibration modes are calculated with the finite displacement method.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/vibration.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.VibrationFeature.html#matlantis_features.features.common.vibration.VibrationFeature.__init__)    def __init__(
            self,
            delta: float = 0.01,
            nfree: int = 2,
            method: str = "standard",  # standard or frederiksen
            direction: str = "central",  # central, forward or backward
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            estimator_fn: Optional[EstimatorFnType] = None,
            num_threads: int = 20,
        ) -> None:
            super(VibrationFeature, self).__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            method = method.lower()
            assert nfree in [2, 4]
            assert method in ["standard", "frederiksen"]
            assert direction in ["central", "forward", "backward"]
            with self.init_scope():
                self.delta = delta
                self.nfree = nfree
                self.method = method
                self.direction = direction
                self.show_progress_bar = show_progress_bar
                self.tqdm_options = tqdm_options
                self.show_logger = show_logger
                self.estimator_fn = estimator_fn
                self.num_threads = num_threads
    
    
    
        __init__.__doc__ = _VIBRATION_DOC_STRING
    
        def __call__(self, atoms: Union[MatlantisAtoms, ASEAtoms], indices: Optional[List[int]] = None) -> VibrationFeatureResult:
            """Run the vibration calculation.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    If 'None' is provided, all atoms are allowed to vibrate. Defaults to None.
            Returns:
                VibrationFeatureResult :
                    The vibration calculation results, which include the force constant, vibration
                    frequencies and vibration modes.
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            vib = FastVibrations(
                delta=self.delta,
                nfree=self.nfree,
                method=self.method,
                direction=self.direction,
                show_progress_bar=self.show_progress_bar,
                tqdm_options=self.tqdm_options,
                show_logger=self.show_logger,
                estimator_fn=self.estimator_fn,
                num_threads=self.num_threads,
            )
            res = vib.run(
                ase_atoms,
                indices=indices,
            )
            with get_calculator(estimator_fn=self.estimator_fn) as calculator:
                ase_atoms.calc = calculator
                potential_energy = ase_atoms.get_potential_energy()
    
            results = VibrationFeatureResult(
                atoms=MatlantisAtoms(ase_atoms.copy()),
                indices=indices,
                delta=self.delta,
                nfree=self.nfree,
                method=self.method,
                direction=self.direction,
                all_forces=res.all_forces,
                force_constant=res.H,
                frequencies=res.frequencies.real,
                image_frequencies=res.frequencies.imag,
                mode=res.mode,
                zero_point_energy=res.zero_point_energy,
                potential_energy=potential_energy,
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
            n_pfp_run = n_atoms * 3 * self.nfree + 1
            return FeatureCost(n_pfp_run=n_pfp_run, time_additional=0.1)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.FastVibrationsResult.html#matlantis_features.features.common.vibration.FastVibrationsResult)@dataclass
    class FastVibrationsResult:
        """A dataclass for result of FastVibration."""
    
        H: np.ndarray
        all_forces: np.ndarray
        frequencies: np.ndarray
        mode: np.ndarray
        zero_point_energy: float
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.FastVibrations.html#matlantis_features.features.common.vibration.FastVibrations)class FastVibrations:
        """The class for calculating the force constant, vibration frequency and vibration mode \
            using the finite displacement method."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.FastVibrations.html#matlantis_features.features.common.vibration.FastVibrations.__init__)    def __init__(
            self,
            delta: float = 0.01,
            nfree: int = 2,
            method: str = "standard",
            direction: str = "central",
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            estimator_fn: Optional[EstimatorFnType] = None,
            num_threads: int = 10,
        ) -> None:
            self.hessian = _Hessian(
                delta,
                nfree,
                method,
                direction,
                show_progress_bar,
                tqdm_options,
                show_logger,
                estimator_fn,
                num_threads,
            )
    
    
    
        __init__.__doc__ = _VIBRATION_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.vibration.FastVibrations.html#matlantis_features.features.common.vibration.FastVibrations.run)    def run(
            self,
            atoms: Union[MatlantisAtoms, ASEAtoms],
            indices: Optional[List[int]] = None,
        ) -> FastVibrationsResult:
            """Run vibration calculation.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    Defaults to None.
            Returns:
                FastVibrationsResult : The vibration calculation result.
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            if indices and (max(indices) >= len(ase_atoms) or min(indices) < 0):
                raise IndexError
    
            if indices is None:
                indices = list(range(len(ase_atoms)))
    
            masses = ase_atoms.get_masses()[indices] if indices else ase_atoms.get_masses()
            displacements = self.hessian._displacements(ase_atoms, indices)
            eq_forces, all_forces = self.hessian._calculate_forces(displacements, indices)
            H = self.hessian._calculate_hessian(eq_forces, all_forces, indices)
            hnu, mode = self._calculate_frequencies(H, masses)
            zpe = float(0.5 * hnu.real.sum())
            res = FastVibrationsResult(
                H=H,
                all_forces=all_forces,
                frequencies=hnu,
                mode=mode,
                zero_point_energy=zpe,
            )
            return res
    
    
    
        def _calculate_frequencies(self, H: np.ndarray, masses: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            """Calculate the vibration frequencies and vibration modes.
    
            Args:
                H (np.ndarray): The force constant.
                masses (np.ndarray): The atomic mass of each atom.
            Returns:
                tuple[np.ndarray, np.ndarray] : The vibration frequencies and vibration modes
            """
            assert len(H) == len(masses) * 3
            if 0 in masses:
                raise RuntimeError(
                    "Zero mass encountered in one or more of the vibrated atoms. Use Atoms.set_masses() to set all masses to non-zero values."
                )
    
            im = np.repeat(masses**-0.5, 3)
            omega2, modes = np.linalg.eigh(im[:, None] * H * im)
            modes = modes.T.copy()
    
            # Conversion factor:
            s = units._hbar * 1e10 / sqrt(units._e * units._amu)
            hnu = s * omega2.astype(complex) ** 0.5
            return hnu, modes
    
    
    
