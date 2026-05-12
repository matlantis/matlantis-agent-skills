# Source code for matlantis_features.features.common.hessian
    
    
    import logging
    import warnings
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional, Tuple, Union
    
    import numpy as np
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    from matlantis_features.utils.parallel import _parallel_calculation
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.hessian.HessianFeatureResult.html#matlantis_features.features.common.hessian.HessianFeatureResult)@dataclass
    class HessianFeatureResult:
        """A dataclass for result of HessianFeature."""
    
        atoms: MatlantisAtoms
        indices: Optional[List[int]]
        delta: float
        nfree: int
        method: str
        direction: str
        H: np.ndarray
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the HessianFeatureResult.
    
            Returns:
                dict[str, Any] : The vibration calculation results in dictionary format.
            """
            return asdict(self)
    
    
    
    
    _HESSIAN_DOC_STRING = """Initialize an instance.
    
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
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.hessian.HessianFeature.html#matlantis_features.features.common.hessian.HessianFeature)class HessianFeature(FeatureBase):
        """The matlantis-feature for calculating molecular vibration.
    
        The vibration properties, including the force constant, vibration frequencies and \
        vibration modes are calculated with the finite displacement method.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.hessian.HessianFeature.html#matlantis_features.features.common.hessian.HessianFeature.__init__)    def __init__(
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
            super(HessianFeature, self).__init__()
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
    
    
    
        __init__.__doc__ = _HESSIAN_DOC_STRING
    
        def __call__(self, atoms: Union[MatlantisAtoms, ASEAtoms], indices: Optional[List[int]] = None) -> HessianFeatureResult:
            """Calculate the hessian matrix.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    If 'None' is provided, all atoms are allowed to vibrate. Defaults to None.
            Returns:
                HessianFeatureResult: The hessian matrix.
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            hessian = _Hessian(
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
            res = hessian.run(
                ase_atoms,
                indices=indices,
            )
    
            results = HessianFeatureResult(
                atoms=MatlantisAtoms(ase_atoms.copy()),
                indices=indices,
                delta=self.delta,
                nfree=self.nfree,
                method=self.method,
                direction=self.direction,
                H=res,
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
    
    
    
    
    class _Hessian:
        """The class for calculating the hessian matrix using the finite displacement method."""
    
        def __init__(
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
            assert nfree in [2, 4]
            assert method in ["standard", "frederiksen"]
            assert direction in ["central", "forward", "backward"]
    
            self.delta = delta
            self.nfree = nfree
            self.method = method
            self.direction = direction
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.estimator_fn = estimator_fn
            self.num_threads = num_threads
    
            if self.nfree == 2:
                self.ndis = [-1, 1]
            else:
                self.ndis = [-2, -1, 1, 2]
    
        __init__.__doc__ = _HESSIAN_DOC_STRING
    
        def run(
            self,
            atoms: Union[MatlantisAtoms, ASEAtoms],
            indices: Optional[List[int]] = None,
        ) -> np.ndarray:
            """Get hessian matrix.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    Defaults to None.
            Returns:
                np.ndarray : The hessian matrix of (3*n_atoms, 3*n_atoms)
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            if len(ase_atoms) > 2000 and self.num_threads > 20:
                warnings.warn(
                    "Internal Server Error may be caused when using many threads (>20) "
                    "for the parallel calculation of the large input structure (>2000 atoms).",
                    UserWarning,
                )
    
            if indices and (max(indices) >= len(ase_atoms) or min(indices) < 0):
                raise IndexError
    
            if indices is None:
                indices = list(range(len(ase_atoms)))
    
            displacements = self._displacements(ase_atoms, indices)
            eq_forces, all_forces = self._calculate_forces(displacements, indices)
            H = self._calculate_hessian(eq_forces, all_forces, indices)
            return H
    
        def _displacements(self, atoms: ASEAtoms, indices: Optional[List[int]] = None) -> List[ASEAtoms]:
            """Create a set of structures with displacements.
    
            Args:
                atoms (ASEAtoms): The input structure.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    Defaults to None.
            Returns:
                list[ASEAtoms] : A set of structures with displacments.
            """
            if indices is None:
                indices = list(range(len(atoms)))
    
            displacements = [atoms.copy()]
            for a in indices:
                for i in range(3):
                    for n in self.ndis:
                        disp_atoms = atoms.copy()
                        disp_atoms.positions[a, i] += n * self.delta
                        displacements.append(disp_atoms)
            return displacements
    
        def _calculate_forces(
            self,
            displacements: List[ASEAtoms],
            indices: Optional[List[int]] = None,
        ) -> Tuple[np.ndarray, np.ndarray]:
            """Calculate the force for the displaced structure in parallel.
    
            Args:
                displacements (list[ASEAtoms]): a set of displaced structures created with the
                    self._displacements method.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    Defaults to None.
            Returns:
                tuple[np.ndarray, np.ndarray] :
                    The forces of the equilibrium structure and displaced structures.
            """
            logger = logging.getLogger(__name__)
    
            if indices is None:
                indices = list(range(len(displacements[0])))
    
            forces: List[np.ndarray] = []
            ind = [-1] + [i for i in indices for _ in range(3 * self.nfree)]
            disp = [np.array([0.0, 0.0, 0.0])] + [
                np.array([self.delta * n if j == i else 0.0 for j in range(3)]) for a in indices for i in range(3) for n in self.ndis
            ]
            results_list = _parallel_calculation(
                displacements,
                estimator_fn=self.estimator_fn,
                num_threads=self.num_threads,
                properties=["energy", "forces"],
                show_progress_bar=self.show_progress_bar,
                tqdm_options=self.tqdm_options,
            )
            for a, d, results in zip(ind, disp, results_list):
                if self.show_logger:
                    if a == -1:
                        logger.info("Force calculation (no displacement)")
                    else:
                        logger.info(f"Force calculation ({a}th atoms displacement by [{d[0]: 5.3f} {d[1]: 5.3f} {d[2]: 5.3f}])")
                f = results["forces"]
                if self.method == "frederiksen" and a > 0:
                    f[a] -= f.sum(0)
                forces.append(f[indices])
    
            eq_forces = forces[0]
            all_forces = np.array(forces[1:]).reshape([len(indices), 3, self.nfree, len(indices), 3])
            return eq_forces, all_forces
    
        def _calculate_hessian(
            self,
            eq_forces: np.ndarray,
            all_forces: np.ndarray,
            indices: Optional[List[int]] = None,
        ) -> np.ndarray:
            """Calculate the force constant from the forces of displaced structures.
    
            Args:
                eq_forces (np.ndarray): The forces of the equilibrium structure.
                all_forces (np.ndarray): The forces of the displaced structures.
                indices (list[int] or None, optional): List of indices of atoms allowed to vibrate.
                    Defaults to None.
            Returns:
                np.ndarray : The force constant.
            """
            if indices is None:
                indices = list(range(len(eq_forces)))
    
            n = 3 * len(indices)
            feq = eq_forces.reshape([1, 1, len(indices), 3])
            H: np.ndarray
            if self.direction == "central":
                if self.nfree == 2:
                    coef = np.array([0.5, -0.5]).reshape([1, 1, 2, 1, 1])
                else:
                    coef = np.array([-1 / 12.0, 8 / 12.0, -8 / 12.0, 1 / 12.0]).reshape([1, 1, 4, 1, 1])
                H = np.sum(all_forces * coef, axis=2).reshape(n, n)
            elif self.direction == "forward":
                H = (feq - all_forces[:, :, self.nfree // 2, :, :]).reshape(n, n)
            else:
                H = (all_forces[:, :, self.nfree // 2 - 1, :, :] - feq).reshape(n, n)
    
            H /= 2 * self.delta
            H += H.copy().T
            return H
    
