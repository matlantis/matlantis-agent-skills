# Source code for matlantis_features.features.phonon.symmetry
    
    
    from typing import List, Optional, Tuple, Union
    
    import numpy as np
    from ase import Atoms
    from spglib import get_symmetry_dataset
    
    # directions = np.eye(3, dtype=int)
    directions = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, -1, 0],
            [1, 0, -1],
            [0, 1, -1],
            [1, 1, 1],
            [1, 1, -1],
            [1, -1, 1],
            [-1, 1, 1],
        ]
    )
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations)class SymmetryOperations:
        """All the symmetry operations in the input structure."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.__init__)    def __init__(self, atoms: Atoms, symprec: float = 1e-5) -> None:
            """Initialize an instance.
    
            Args:
                atoms (Atoms): The input structure.
                symprec (float, optional): Distance tolerance to find crystal symmetry. Defaults to 1e-05.
            """
            spglib_cell = (
                atoms.get_cell()[:],
                atoms.get_scaled_positions(),
                atoms.get_atomic_numbers(),
            )
            self.atoms = atoms
            self.positions = self.atoms.get_scaled_positions()
            self.cell = atoms.get_cell()[:]
            self.symprec = symprec
            self.database: np.ndarray = get_symmetry_dataset(spglib_cell, symprec)
            self.num_symmetry_operations = len(self.database["rotations"])
            self._permutations: List[Optional[np.ndarray]] = [None for _ in range(self.num_symmetry_operations)]
            self._rotation_coord: List[Optional[np.ndarray]] = [None for _ in range(self.num_symmetry_operations)]
            self._rotation_on_atom: List[Optional[Tuple[np.ndarray, np.ndarray]]] = [None for _ in range(len(self.atoms))]
    
    
    
        @property
        def rotations(self) -> np.ndarray:
            """Get the rotation matrix of all the symmetry operations.
    
            Returns:
                np.ndarray : rotation matrix.
            """
            rotations: np.ndarray = self.database["rotations"]
            return rotations
    
        @property
        def translations(self) -> np.ndarray:
            """Get the translation vector of all the symmetry operations.
    
            Returns:
                np.ndarray : translation vector.
            """
            translations: np.ndarray = self.database["translations"]
            return translations
    
        @property
        def equivalent_atoms(self) -> np.ndarray:
            """Mapping the atoms to the independent atoms.
    
            Returns:
                np.ndarray : mapping index.
            """
            equivalent_atoms: np.ndarray = self.database["equivalent_atoms"]
            return equivalent_atoms
    
        @property
        def independent_atoms(self) -> np.ndarray:
            """List of independent atoms.
    
            Returns:
                np.ndarray : List of independent atoms.
            """
            independent_atoms: np.ndarray = np.unique(self.database["equivalent_atoms"])
            return independent_atoms
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.get_rotations_coord)    def get_rotations_coord(self, index: int) -> np.ndarray:
            """The rotation matrix for the cartesian coordinations instead of fractional coordinations.
    
            Args:
                index (int): The index of rotation operation.
            Returns:
                np.ndarray : The rotation matrix for the cartesian coordinations.
            """
            if self._rotation_coord[index] is None:
                r = self.database["rotations"][index]
                self._rotation_coord[index] = similarity_transformation(self.cell.T, r)
            return self._rotation_coord[index]  # type: ignore
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.get_rotation_on_atom)    def get_rotation_on_atom(self, index: int, return_indices: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
            """Find the rotation matrix center on the given atom.
    
            Args:
                index (int): The index of the center atom.
                return_indices (bool, optional): If True, also return the indices of rotation matrix.
                    Defaults to False.
            Returns:
                np.ndarray or tuple[np.ndarray, np.ndarray] : The rotation matrix center on the given atom.
            """
            if self._rotation_on_atom[index] is None:
                pos = self.positions[index]
                rotation_on_atom = []
                indices = []
                for i, (r, t) in enumerate(zip(self.database["rotations"], self.database["translations"])):
                    rotated_pos = np.dot(pos, r.T) + t
                    diff = rotated_pos - pos
                    diff -= np.rint(diff)
                    diff = np.dot(diff, self.cell)
                    if np.sum(diff**2) < self.symprec**2:
                        rotation_on_atom.append(r)
                        indices.append(i)
                self._rotation_on_atom[index] = (np.array(rotation_on_atom), np.array(indices))
            if return_indices:
                return self._rotation_on_atom[index]  # type: ignore
            else:
                return self._rotation_on_atom[index][0]  # type: ignore
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.get_permutations_on_atom)    def get_permutations_on_atom(self, index: int) -> np.ndarray:
            """If the initial strucutre R = {R_i} was rotated with matrix {r} on the atoms[index] and\
            get the rotated structure R' = {R_i}', we can get a serial of index from this function\
            that makes R'[ind] = R.
    
            Args:
                index (int): the index of the center atom.
            Returns:
                np.ndarray : the permutations index.
            """
            rotations = self.get_rotation_on_atom(index)
            positions = self.positions - self.positions[index]
            return np.array(  # type: ignore
                [get_permutations(positions.dot(r.T), positions) for r in rotations]
            )
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.get_translations_on_atom)    def get_translations_on_atom(self, index: int) -> int:
            """To find a set of symmetry operation (both rotation and translation) that does not change \
            the position of atoms[index].
    
            Args:
                index (int): the index of the unmoved atom
            Returns:
                int : the index of rotation and translation operations.
            """
            for i in range(self.num_symmetry_operations):
                p = self.get_permutations(i)
                if p[index] in self.independent_atoms:
                    return i
            raise ValueError
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations.get_permutations)    def get_permutations(self, index: int) -> np.ndarray:
            """If the initial structure R is transformed with rotation[index] and translation[index]\
            and get the rotated structure R'. we can get a serial of index from this function which\
            enable R' = R[ind].
    
            Args:
                index (int): the index of symmetry operation.
            Returns:
                np.ndarray : the permutations index.
            """
            if self._permutations[index] is None:
                positions = self.atoms.get_scaled_positions()
                r = self.database["rotations"][index]
                t = self.database["translations"][index]
                rotated_positions = np.dot(positions, r.T) + t
                rotated_positions -= np.floor(rotated_positions)
                p = get_permutations(positions, rotated_positions)
                self._permutations[index] = p
            return self._permutations[index]  # type: ignore
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry.html#matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry)def get_least_displacements_from_symmetry(
        atoms: Atoms, atoms_N: Atoms, symmetry_operations: SymmetryOperations, delta: float
    ) -> Tuple[List[Atoms], np.ndarray, np.ndarray]:
        """To find the least number of displacements that need to calculate force constant.
    
        Args:
            atoms (Atoms): The input structure.
            atoms_N (Atoms): The supercell.
            symmetry_operations (SymmetryOperations): The symmetry operations in the supercell.
            delta (float): The distance of finite displacement applied to each atom.
        Returns:
            tuple[list[Atoms], np.ndarray, np.ndarray] :
                A list of displaced structures, the index of the displaced atom, and the displacment
                vector of the displaced atom.
        """
        disp_atoms: List[Atoms] = []
        disp_index: List[int] = []
        disp: List[np.ndarray] = []
        for i in symmetry_operations.independent_atoms:
            direc_plus: List[np.ndarray] = []
            direc_minus: List[np.ndarray] = []
            rotation_on_atom: np.ndarray = symmetry_operations.get_rotation_on_atom(i)  # type: ignore
            found = _get_directions_1(direc_plus, rotation_on_atom)
            if not found:
                found = _get_directions_2(direc_plus, rotation_on_atom)
            if not found:
                direc_plus += [directions[0], directions[1], directions[2]]
    
            for d in direc_plus:
                found_minus = _get_minus_direction(d, rotation_on_atom)
                if not found_minus:
                    direc_minus.append(-d)
    
            direc = direc_plus + direc_minus
            for d in direc:
                atoms_N_copy = atoms_N.copy()
                atoms_N_copy.positions[i] += d * delta / np.linalg.norm(d)
                disp_atoms.append(atoms_N_copy)
                disp_index.append(i)
                disp.append(d * delta / np.linalg.norm(d))
        return disp_atoms, np.array(disp_index), np.array(disp)
    
    
    
    
    def _get_directions_1(direc_plus: List[np.ndarray], rotation_on_atom: np.ndarray) -> bool:
        num_sym = len(rotation_on_atom)
        for d in directions:
            for j in range(num_sym):
                for k in range(j + 1, num_sym):
                    det = np.linalg.det(np.array([d, np.dot(d, rotation_on_atom[j].T), np.dot(d, rotation_on_atom[k].T)]))
                    if det != 0:
                        direc_plus.append(d)
                        return True
        return False
    
    
    def _get_directions_2(direc_plus: List[np.ndarray], rotation_on_atom: np.ndarray) -> bool:
        num_sym = len(rotation_on_atom)
        for d in directions:
            for j in range(num_sym):
                for d2 in directions:
                    det = np.linalg.det(np.array([d, np.dot(d, rotation_on_atom[j].T), d2]))
                    if det != 0:
                        direc_plus += [d, d2]
                        return True
        return False
    
    
    def _get_minus_direction(d: np.ndarray, rotation_on_atom: np.ndarray) -> bool:
        num_sym = len(rotation_on_atom)
        for j in range(num_sym):
            rotated_d = -np.dot(d, rotation_on_atom[j].T)
            if np.all(d == rotated_d):
                return True
        return False
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.resolve_force_constants.html#matlantis_features.features.phonon.symmetry.resolve_force_constants)def resolve_force_constants(
        atoms: Atoms,
        atoms_N: Atoms,
        symmetry_operations: SymmetryOperations,
        forces: np.ndarray,
        disp_index: np.ndarray,
        disp: np.ndarray,
    ) -> np.ndarray:
        """Calculate the force constant from the least number of displacemented structures.
    
        Args:
            atoms (Atoms): The input structure.
            atoms_N (Atoms): The supercell.
            symmetry_operations (SymmetryOperations): The symmetry operations in the supercell.
            forces (np.ndarray): The forces of the displaced structures.
            disp_index (np.ndarray): The index of the displaced atom.
            disp (np.ndarray): the displacment vector of the displaced atom.
        Returns:
            np.ndarray : The force constant.
        """
        force_constants = np.zeros((len(atoms), len(atoms_N), 3, 3), dtype="double", order="C")
        disp_index_unique = np.unique(disp_index)
        for i in disp_index_unique:
            disps_on_atom = disp[disp_index == i]
            forces_i = forces[disp_index == i]
            rotation_on_atom, operation_indices = symmetry_operations.get_rotation_on_atom(i, return_indices=True)
            rotation_coords = np.array([symmetry_operations.get_rotations_coord(i) for i in operation_indices])
            permutation_on_atom = symmetry_operations.get_permutations_on_atom(i)
            rotated_disps = np.array([np.dot(r, u) for u in disps_on_atom for r in rotation_coords])
            inv_displacements = np.linalg.pinv(rotated_disps)
            for j in range(len(atoms_N)):
                rotated_forces_list = []
                for f in forces_i:
                    f_perm = f[permutation_on_atom[:, j]]
                    for fp, r in zip(f_perm, rotation_coords):
                        rotated_forces_list.append(np.dot(fp, r.T))
                rotated_forces = np.array(rotated_forces_list)
                force_constants[i][j] = -np.dot(inv_displacements, rotated_forces)
        map_atoms = symmetry_operations.equivalent_atoms
        for i in range(len(atoms)):
            if i not in symmetry_operations.independent_atoms:
                map_syms = symmetry_operations.get_translations_on_atom(i)
                r_cart = symmetry_operations.get_rotations_coord(map_syms)
                permutation = symmetry_operations.get_permutations(map_syms)
                for j in range(len(atoms_N)):
                    fc_computed = force_constants[map_atoms[i]][permutation[j]]
                    force_constants[i][j] += np.dot(r_cart.T, np.dot(fc_computed, r_cart))
        return force_constants  # type: ignore
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.similarity_transformation.html#matlantis_features.features.phonon.symmetry.similarity_transformation)def similarity_transformation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Calculate the similarity matrix, B = P^-1 A P.
    
        Args:
            a (np.ndarray): the change of basis matrix.
            b (np.ndarray): the input matrix.
        Returns:
            np.ndarray : Transformed matrix.
        """
        return np.dot(a, np.dot(b, np.linalg.inv(a)))  # type: ignore
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.symmetry.get_permutations.html#matlantis_features.features.phonon.symmetry.get_permutations)def get_permutations(pos_a: np.ndarray, pos_b: np.ndarray, symprec: float = 1e-5) -> np.ndarray:
        """Get the permuation index which enables pos_b = pos_a[ind].
    
        Args:
            pos_a (np.ndarray): coordinations before permutation.
            pos_b (np.ndarray): coordinations after permutation.
            symprec (float, optional): distance tolerance to find crystal symmetry. Defaults to 1e-05.
        Returns:
            np.ndarray : the permuation index.
        """
        diff = pos_a[np.newaxis] - pos_b[:, np.newaxis]
        diff -= np.rint(diff)
        d = np.sum(diff**2, axis=2)
        row, ind = np.where(d < symprec**2)
        if not np.all(row == np.arange(len(pos_a))) or len(ind) != len(pos_a):
            print(row)
            raise PermutationIndexError
        return ind  # type: ignore
    
    
    
    
    class PermutationIndexError(Exception):
        """Raise error when cannot found the permutation index."""
    
        pass
    
