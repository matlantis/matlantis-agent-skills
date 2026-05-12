# Source code for matlantis_features.features.phonon.force_constant
    
    
    import logging
    from dataclasses import dataclass
    from typing import Any, Dict, List, Optional, Tuple, Union
    
    import numpy as np
    from ase import Atoms
    from ase.calculators.calculator import Calculator
    from tqdm.auto import tqdm
    from tqdm.contrib.logging import logging_redirect_tqdm
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.symmetry import (
        SymmetryOperations,
        get_least_displacements_from_symmetry,
        resolve_force_constants,
    )
    from matlantis_features.features.phonon.utils import get_primitive_structure
    from matlantis_features.utils import FeatureCost, get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult)@dataclass
    class ForceConstantFeatureResult:
        """A dataclass for result of ForceConstantFeature."""
    
        force_constant: np.ndarray
        supercell: np.ndarray
        delta: float
        name: str
        unit_cell_atoms: MatlantisAtoms
        unit_cell_potential_energy: float
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.ForceConstantFeature.html#matlantis_features.features.phonon.force_constant.ForceConstantFeature)class ForceConstantFeature(FeatureBase):
        """The matlantis-feature for calculating the force constant for the given structure.
    
        The force constant, which is the basis of phonon analysis, will be calculated with the \
        finite displacement method, in which the second derivative of the potential energy is \
        estimated by displacing each atom with a small distance within an enlarged supercell.
        The results will be used in the post analysis, such as PostPhononBandFeature.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/force_constant.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.ForceConstantFeature.html#matlantis_features.features.phonon.force_constant.ForceConstantFeature.__init__)    def __init__(
            self,
            supercell: Tuple[int, int, int],
            delta: float,
            method: str = "Frederiksen",
            direction: str = "central",
            symmetrize: int = 3,
            acoustic: bool = True,
            cutoff: Optional[float] = None,
            use_symmetry: bool = False,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            estimator_fn: Optional[EstimatorFnType] = None,
        ):
            """Initialize an instance.
    
            Args:
                supercell (tuple[int, int, int]): The enlarged supercell that is created by elongating
                    the input structure by (n, m, l) times along the axes.
                delta (float): The distance of finite displacement applied to each atom. The unit
                    is Angstrom.
                method (str, optional): The method to calculate force constant from the displaced
                    supercells and their forces. This parameter can be 'Standard' or 'Frederiksen'.
                    Defaults to 'Frederiksen'.
                direction (str, optional): This parameter discribe how the displacment was found. Must
                    be in 'forward', 'backward' or 'central'. When 'direction=forward'/
                    'direction=backward', each atom is displaced along the positive/negative direction
                    of the axes. When 'direction=center', each atom is displaced along both positive
                    and negative directions of the axes. Defaults to 'central'.
                symmetrize (int, optional): After restoring the acoustic sum rule, the symmetry of the force
                    constant is broken. The parameter specifies the number of interactions to symmetrize
                    the force constant. Defaults to 3.
                acoustic (bool, optional): Restores the acoustic sum rule of the force constant. Defaults
                    to True.
                cutoff (float or None, optional): The corresponding elements in the force constant are set
                    to 0 if the interatomic distance between atoms is larger than cutoff. If None is
                    provided, the full force constant will be used. Defaults to None.
                use_symmetry (bool): If True, the symmetry of input structure will be used to reduce the
                    number of force estimations. Defaults to False.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            super(ForceConstantFeature, self).__init__()
            self.check_estimator_fn(estimator_fn)
            method = method.lower()
            direction = direction.lower()
    
            assert method in ["frederiksen", "standard"]
            assert direction in ["central", "forward", "backward"]
    
            with self.init_scope():
                self.supercell = supercell
                self.delta = delta
                self.method = method
                self.direction = direction
                self.symmetrize = symmetrize
                self.acoustic = acoustic
                self.cutoff = cutoff
                self.use_symmetry = use_symmetry
                self.show_progress_bar = show_progress_bar
                self.tqdm_options = tqdm_options
                self.show_logger = show_logger
                self.estimator_fn = estimator_fn
    
    
    
        def __call__(
            self,
            atoms: Union[MatlantisAtoms, ASEAtoms],
            primitive_matrix: Optional[Union[np.ndarray, str]] = None,
            prec: float = 1e-3,
        ) -> ForceConstantFeatureResult:
            """Calculates force constant for the input structure.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
                primitive_matrix (np.ndarray or str or None, optional):
                    When 'primitive_matrix' is not None,
                    the input structure will be transformed to the primitive cell according to the
                    axes provided here. Then, the following force constant calculation will be
                    performed with the reduced primitive cell. Defaults to None.
                prec (float, optional): Wrapper the input structure into the primitive cell with certain
                    precision. The distance between the wrapped atom and the atom in the primitive cell
                    should be smaller than prec. Defaults to 0.001.
            Returns:
                ForceConstantFeatureResult: The force constant.
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            if primitive_matrix is not None:
                primitive_atoms = get_primitive_structure(ase_atoms, primitive_matrix, prec)
            else:
                primitive_atoms = ase_atoms
    
            with get_calculator(estimator_fn=self.estimator_fn) as calculator:
                primitive_atoms.calc = calculator
    
                phonons = FastPhonons(
                    calculator,
                    self.supercell,
                    self.delta,
                    self.method,
                    self.direction,
                    self.symmetrize,
                    self.acoustic,
                    self.use_symmetry,
                    self.cutoff,
                    self.show_progress_bar,
                    self.tqdm_options,
                    self.show_logger,
                )
                res = phonons.run(primitive_atoms)
                unit_cell_potential_energy = primitive_atoms.get_potential_energy()
    
            return ForceConstantFeatureResult(
                force_constant=res,
                supercell=np.array(self.supercell),
                delta=self.delta,
                name=self.name,
                unit_cell_atoms=MatlantisAtoms(primitive_atoms.copy()),
                unit_cell_potential_energy=unit_cell_potential_energy,
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): The number of atoms.
            Returns:
                FeatureCost :
                    A dataclass containing cost information for the feature.
            """
            n_pfn_run = n_atoms * 6 + 1
            return FeatureCost(n_pfp_run=n_pfn_run, time_additional=1.0)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.FastPhonons.html#matlantis_features.features.phonon.force_constant.FastPhonons)class FastPhonons(object):
        """Class for calculating force constant using the finite displacement method."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.FastPhonons.html#matlantis_features.features.phonon.force_constant.FastPhonons.__init__)    def __init__(
            self,
            calculator: Calculator,
            supercell: Tuple[int, int, int],
            delta: float = 0.1,
            method: str = "frederiksen",
            direction: str = "central",
            symmetrize: int = 3,
            acoustic: bool = True,
            use_symmetry: bool = False,
            cutoff: Optional[float] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
        ):
            """Initialize an instance.
    
            Args:
                calculator (Calculator): The calculator to estimate energy and force.
                supercell (tuple[int, int, int]): The number of repetitions of input structure.
                delta (float, optional): The distance of finite displacment. Defaults to 0.1.
                method (str, optional): The method to calculate force constant. Defaults to
                    'frederiksen'.
                direction (str, optional): How the displacement is found. Defaults to 'central'.
                symmetrize (int, optional): The number of interactions to symmetrize the force
                    constant. Defaults to 3.
                acoustic (bool, optional): Restores the acoustic sum rule of the force constant.
                    Defaults to True.
                use_symmetry (bool): If True, the symmetry of input structure will be used to reduce the
                    number of force estimations. Defaults to False.
                cutoff (float or None, optional):  The corresponding elements in the force constant is
                    set to 0 if the interatomic distance between atoms is larger than cutoff. Defaults
                    to None.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
            """
            method = method.lower()
            direction = direction.lower()
            assert method in ["standard", "frederiksen"]
            assert direction in ["central", "forward", "backward"]
    
            self.calc = calculator
            self.supercell = supercell
            self.delta = delta
            self.method = method
            self.direction = direction
            self.symmetrize = symmetrize
            self.acoustic = acoustic
            self.use_symmetry = use_symmetry
            self.cutoff = cutoff
            self.n_supercell = np.prod(supercell)
            self.n_free = 2 if direction == "central" else 1
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.force_constant.FastPhonons.html#matlantis_features.features.phonon.force_constant.FastPhonons.run)    def run(
            self,
            atoms: Union[MatlantisAtoms, ASEAtoms],
        ) -> np.ndarray:
            """Runs the force constant calculation.
    
            Args:
                atoms (MatlantisAtoms or ASEAtoms): The input structure.
            Returns:
                np.ndarray : The force constant
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
    
            atoms_N = ase_atoms * self.supercell
            symmetry_operations: Optional[SymmetryOperations]
            if self.use_symmetry:
                symmetry_operations = SymmetryOperations(atoms_N, symprec=1e-5)
            else:
                symmetry_operations = None
    
            disp_atoms, disp_index, disp = self._displacements(ase_atoms, atoms_N, symmetry_operations)
            forces = self._calculate_forces(ase_atoms, atoms_N, disp_atoms, disp_index, disp)
            force_constant = self._calculate_force_constant(ase_atoms, atoms_N, symmetry_operations, forces, disp_index, disp)
            return force_constant
    
    
    
        def _displacements(
            self, atoms: ASEAtoms, atoms_N: ASEAtoms, symmetry_operations: Optional[SymmetryOperations]
        ) -> Tuple[List[ASEAtoms], np.ndarray, np.ndarray]:
            """Creates a set of supercells with displacements.
    
            Args:
                atoms (ASEAtoms): The input primitive cell structure.
                atoms_N (ASEAtoms): The supercell structure.
                symmetry_operations (Optional[SymmetryOperations]): The symmetry operations in the supercell.
            Returns:
                tuple[list[ASEAtoms], np.ndarray, np.ndarray] :
                    A list of displaced structures, the index of the displaced atom, and the displacment
                    vector of the displaced atom.
            """
            if not self.use_symmetry:
                if self.direction == "central":
                    ndis = [-1, 1]
                elif self.direction == "forward":
                    ndis = [1]
                else:
                    ndis = [-1]
    
                disp_atoms = [atoms_N.copy()]
                disp_index_list = [-1]
                disp_direction_list = [np.array([0, 0, 0])]
                for a in range(len(atoms)):
                    for i in range(3):
                        direc = np.zeros(3, dtype=int)
                        direc[i] = 1
                        for n in ndis:
                            disp_atoms_N = atoms_N.copy()
                            disp_atoms_N.positions[a] += n * direc * self.delta
                            disp_atoms.append(disp_atoms_N)
                            disp_index_list.append(a)
                            disp_direction_list.append(n * direc * self.delta)
                return disp_atoms, np.array(disp_index_list), np.array(disp_direction_list)
            else:
                disp_atoms, disp_index, disp_direction = get_least_displacements_from_symmetry(
                    atoms,
                    atoms_N,
                    symmetry_operations,  # type: ignore
                    self.delta,  # type: ignore
                )
    
            return disp_atoms, disp_index, disp_direction
    
        def _calculate_forces(
            self,
            atoms: ASEAtoms,
            atoms_N: ASEAtoms,
            disp_atoms: List[ASEAtoms],
            disp_index: np.ndarray,
            disp_direction: np.ndarray,
        ) -> np.ndarray:
            """Calculate the force for the supercells with displacements.
    
            Args:
                atoms (ASEAtoms): The input primitive cell structure.
                atoms_N (ASEAtoms): The supercell structure.
                disp_atoms (list[ASEAtoms]): A set of supercells with displacements created with
                    self._displacements method.
                disp_index (np.ndarray): The index of displaced atom.
                disp_direction (np.ndarray): The vector of atom displacement.
            Returns:
                np.ndarray:
                    The forces of the equilibrium supercell and displaced supercells.
            """
            logger = logging.getLogger(__name__)
    
            pbar = None
            if self.show_progress_bar:
                ntotal = len(disp_atoms)
                if self.tqdm_options is None:
                    pbar = tqdm(total=ntotal)
                else:
                    pbar = tqdm(total=ntotal, **self.tqdm_options)
    
            res_forces: List[np.ndarray] = []
            for atoms_dis, a, disp in zip(disp_atoms, disp_index, disp_direction):
                with logging_redirect_tqdm():
                    if self.show_logger:
                        if a == -1:
                            logger.info("Force calculation (no displacement)")
                        else:
                            logger.info(f"Force calculation ({a}th atoms displacement by [{disp[0]: 5.3f} {disp[1]: 5.3f} {disp[2]: 5.3f}])")
                    if pbar is not None:
                        pbar.update(1)
                f = self.calc.get_forces(atoms_dis)
                res_forces.append(f)
    
            if pbar is not None:
                pbar.close()
            return np.array(res_forces)  # type: ignore
    
        def _calculate_force_constant(
            self,
            atoms: ASEAtoms,
            atoms_N: ASEAtoms,
            symmetry_operations: Optional[SymmetryOperations],
            forces: np.ndarray,
            disp_index: np.ndarray,
            disp: np.ndarray,
        ) -> np.ndarray:
            """Calculates the force constant from the forces of displaced supercells.
    
            Args:
                atoms (ASEAtoms): The input primitive cell structure.
                atoms_N (ASEAtoms): The supercell structure.
                symmetry_operations (Optional[SymmetryOperations]): The symmetry operations in the supercell.
                forces (np.ndarray): The forces of the displaced supercells.
                disp_index (np.ndarray): The index of the displaced atom.
                disp (np.ndarray): the displacment vector of the displaced atom.
            Returns:
                np.ndarray : The force constant.
            """
            if not self.use_symmetry:
                return self._calculate_force_constant_wo_symmetry(atoms, atoms_N, forces)
            else:
                return self._calculate_force_constant_symmetry(
                    atoms,
                    atoms_N,
                    symmetry_operations,  # type: ignore
                    forces,
                    disp_index,
                    disp,  # type: ignore
                )
    
        def _calculate_force_constant_wo_symmetry(self, atoms: Atoms, atoms_N: Atoms, forces: np.ndarray) -> np.ndarray:
            ind = [-1] + [i for i in range(len(atoms)) for _ in range(3 * self.n_free)]
            if self.method == "frederiksen":
                for a, f in zip(ind, forces):
                    if a > 0:
                        f[a] -= f.sum(0)
    
            eq_forces = forces[0]
            all_forces = np.array(forces[1:]).reshape([len(atoms), 3, self.n_free, len(atoms_N), 3])
    
            n = len(atoms) * 3
            feq = eq_forces.reshape([1, 1, len(atoms_N), 3])
            if self.direction == "central":
                coef = np.array([0.5, -0.5]).reshape([1, 1, 2, 1, 1])
                H = np.sum(all_forces * coef, axis=2)
            elif self.direction == "forward":
                H = (feq - all_forces[:, :, 0, :, :]).reshape(n, self.n_supercell, n)
            else:
                H = (all_forces[:, :, 0, :, :] - feq).reshape(n, self.n_supercell, n)
    
            H /= self.delta
            C = np.transpose(H.reshape(n, self.n_supercell, n), [1, 0, 2])
            if self.acoustic:
                if self.symmetrize > 0:
                    # Restore accoustic rule and symmetrize it
                    for _ in range(self.symmetrize):
                        C = _symmetrize_force_constant(C, self.supercell)
                        C = _acoustic(C)
                else:
                    # Restore accoustic rule
                    C = _acoustic(C)
            else:
                if self.symmetrize > 0:
                    # Symmetrize the force constant
                    C = _symmetrize_force_constant(C, self.supercell)
            return C
    
        def _calculate_force_constant_symmetry(
            self,
            atoms: Atoms,
            atoms_N: Atoms,
            symmetry_operations: SymmetryOperations,
            forces: np.ndarray,
            disp_index: np.ndarray,
            disp: np.ndarray,
        ) -> np.ndarray:
            C = resolve_force_constants(atoms, atoms_N, symmetry_operations, forces, disp_index, disp)
            C = np.transpose(C, [0, 2, 1, 3]).reshape(len(atoms) * 3, -1, len(atoms) * 3)
            return np.transpose(C, [1, 0, 2])
    
    
    
    
    def _symmetrize_force_constant(C: np.ndarray, supercell: Tuple[int, int, int]) -> np.ndarray:
        """Symmetrizes the force constant matrix.
    
        Args:
            C (np.ndarray): The force constant.
            supercell (tuple[int, int, int]): The number of repetitions of input structure.
        Returns:
            np.ndarray : Symmetrized force constant.
        """
        x, y = _get_supercell_indices(supercell)
        C = (C[x] + C[y].transpose(0, 2, 1)) / 2.0
        return C
    
    
    def _acoustic(C: np.ndarray) -> np.ndarray:
        """Restores acoustic sum rule of force constants. \
        Let `alpha` and `beta` be one of the 3D coordinates of i-th atom and j-th atom, respectively. \
        For each alpha, beta, and i, it should follow `sum_{L,j}C_{alpha i,beta j}(R_L) = 0` where \
        `R_L` denotes a unit cell.
    
        Args:
            C (np.ndarray): The force constant.
        Returns:
            np.ndarray : The force constant with the restored acoustic sum rule.
        """
        n_atoms = C.shape[1] // 3
        C_reshaped = C.reshape(-1, n_atoms, 3, n_atoms, 3)
        # To be sum to zero among the supercell axis and the axis for `j`
        C_norm: np.ndarray
        C_norm = C_reshaped - np.mean(C_reshaped, axis=(0, 3), keepdims=True)
        return C_norm.reshape(-1, n_atoms * 3, n_atoms * 3)
    
    
    def _get_supercell_indices(s: Tuple[int, int, int]) -> Tuple[np.ndarray, np.ndarray]:
        # Find the pair of unit cells that are symmetric about the origin unit cell.
        s_indices = np.array([[i, j, k] for i in range(s[0]) for j in range(s[1]) for k in range(s[2])])
        ind = s_indices.reshape(-1, 1, 3) + s_indices.reshape(1, -1, 3)
        x, y = np.where(
            ((ind[:, :, 0] == s[0]) | (ind[:, :, 0] == 0))
            & ((ind[:, :, 1] == s[1]) | (ind[:, :, 1] == 0))
            & ((ind[:, :, 2] == s[2]) | (ind[:, :, 2] == 0))
        )
        assert len(x) == len(s_indices)
        return x, y
    
