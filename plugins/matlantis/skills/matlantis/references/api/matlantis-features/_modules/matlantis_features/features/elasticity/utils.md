# Source code for matlantis_features.features.elasticity.utils
    
    
    import math
    import warnings
    from itertools import permutations
    from typing import Callable, Dict, List, Tuple, TypeVar, Union
    
    import numpy as np
    from ase.cell import Cell
    
    from matlantis_features.atoms import ASEAtoms
    
    _T = TypeVar("_T")
    _U = TypeVar("_U")
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.map2d.html#matlantis_features.features.elasticity.utils.map2d)def map2d(data: List[List[_T]], func: Callable[[_T], _U]) -> List[List[_U]]:
        """Apply function to each element in the 2D list.
    
        Args:
            data (list[list[_T]]): Target 2D list.
            func (_T -> _U): function object to apply.
        Returns:
            list[list[_U]] : Resulting 2D list.
        """
        return [[func(i) for i in data_j] for data_j in data]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.select.html#matlantis_features.features.elasticity.utils.select)def select(cond: bool, value_for_true: _T, value_for_false: _T) -> _T:
        """Select 2nd or 3rd arg's value depending on the cond arg.
    
        Args:
            cond (bool): Conditional value.
            value_for_true (_T): Value if cond is True.
            value_for_false (_T): Value if cond is False.
        Returns:
            _T : Result of select operation.
        """
        if cond:
            return value_for_true
        else:
            return value_for_false
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.swap_ij.html#matlantis_features.features.elasticity.utils.swap_ij)def swap_ij(ind: Tuple[_T, _T]) -> Tuple[_T, _T]:
        """Swap the first and second args.
    
        Args:
            ind (tuple[_T, _T]): Input args.
        Returns:
            tuple[_T, _T] : Swapped args.
        """
        return ind[1], ind[0]
    
    
    
    
    voigt_inds = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.voigt_to_ij.html#matlantis_features.features.elasticity.utils.voigt_to_ij)def voigt_to_ij(voigt: int) -> Tuple[int, int]:
        """Convert Voigt notation to 3x3 matrix indices.
    
        Args:
            voigt (int): Voigt index (0-5).
        Returns:
            tuple[int, int] : 2D matrix indices.
        """
        return voigt_inds[voigt]
    
    
    
    
    reverse_voigt_map = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.ij_to_voigt.html#matlantis_features.features.elasticity.utils.ij_to_voigt)def ij_to_voigt(i: int, j: int) -> int:
        """Convert 3x3 matrix indices to Voigt notation.
    
        Args:
            i (int): Matrix index (0-3).
            j (int): Matrix index (0-3).
        Returns:
            int : Voigt index.
        """
        return reverse_voigt_map[i][j]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.is_voigt_diag.html#matlantis_features.features.elasticity.utils.is_voigt_diag)def is_voigt_diag(vix: int) -> bool:
        """Check the Voigt index is diagonal element.
    
        Args:
            vix (int): Voigt index to be checked.
        Returns:
            bool : True if the index indicates the diagonal element.
        """
        if vix >= 6:
            raise IndexError(f"invalid voigt index: {vix}")
        return vix < 3
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.to_unit_vec.html#matlantis_features.features.elasticity.utils.to_unit_vec)def to_unit_vec(a: np.ndarray) -> np.ndarray:
        """Normalize the 3D vector.
    
        Args:
            a (np.ndarray): Input 3D vector of (3,) shape.
        Returns:
            np.ndarray : Normarized 3D vector.
        """
        if len(a) != 3:
            raise ValueError(f"invalid input shape: {a.shape}")
    
        length = np.linalg.norm(a)
        if length < 1e-8:
            raise ValueError(f"The length of input vection is almost zero {length}!")
        return a / length  # type: ignore
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.utils.get_symmetry_lattice.html#matlantis_features.features.elasticity.utils.get_symmetry_lattice)def get_symmetry_lattice(atoms: ASEAtoms) -> ASEAtoms:
        """Rotate the structure to certain oritentation to get the irreducible elastic tensor.
    
        Args:
            atoms (ASEAtoms): The input atoms.
        Returns:
            ASEAtoms : The atoms after rotation.
        """
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        lattice = cell.get_bravais_lattice()
        lattice_system = lattice.lattice_system
        if lattice_system == "triclinic":
            return atoms.copy()
        if lattice_system == "monoclinic":
            typical_cell = _typical_monoclinic_lattice(cell.cellpar())
        else:
            typical_cell = lattice.tocell()
        rotation = _get_rotation_matrix(cell, typical_cell)
        if rotation is False:
            warnings.warn("Cannot get the irreducible elastic tensor.", UserWarning)
            return atoms.copy()
        atoms_new = atoms.copy()
        atoms_new.set_cell(cell[:].dot(rotation))
        atoms_new.set_positions(positions.dot(rotation))
        if lattice_system == "rhombohedral":
            v = atoms_new.cell[:].sum(axis=0)
            atoms_new.rotate(v, [0, 0, 1], rotate_cell=True)
        return atoms_new
    
    
    
    
    def _get_rotation_matrix(cell: Cell, typical_cell: Cell) -> Union[np.ndarray, bool]:
        for perm in permutations(range(3)):
            target_cell = typical_cell[np.array(perm)]
            cell_inv = cell.reciprocal()
            rotation: np.ndarray = cell_inv.T.dot(target_cell[:])
            if _check_rotation_matrix(rotation):
                return rotation
        return False
    
    
    def _check_rotation_matrix(rotation: np.ndarray) -> bool:
        determinant: float = np.linalg.det(rotation)
        det_condition = math.isclose(determinant, 1.0, abs_tol=1e-4) or math.isclose(determinant, -1.0, abs_tol=1e-4)
        return np.allclose(rotation.dot(rotation.T), np.eye(3), atol=1e-4) and det_condition
    
    
    def _typical_monoclinic_lattice(cellpar: Tuple[float, float, float, float, float, float]) -> np.ndarray:
        a, b, c, alpha, beta, gamma = cellpar
        if alpha != 90:
            A, B, C, Alpha = (a, b, c, alpha)
        elif beta != 90:
            A, B, C, Alpha = (b, a, c, beta)
        else:
            A, B, C, Alpha = (c, a, b, gamma)
        return np.array(  # type: ignore
            [
                [A, 0.0, 0.0],
                [0.0, B, 0.0],
                [0.0, np.cos(np.deg2rad(Alpha)) * C, np.sin(np.deg2rad(Alpha)) * C],
            ]
        )
    
    
    zeros_elements: Dict[str, List[List[int]]] = {
        "cubic": [
            [0, 3],
            [0, 4],
            [0, 5],
            [1, 3],
            [1, 4],
            [1, 5],
            [2, 3],
            [2, 4],
            [2, 5],
            [3, 4],
            [3, 5],
            [4, 5],
        ],
        "tetragonal": [[0, 3], [0, 4], [1, 3], [1, 4], [2, 3], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]],
        "orthorhombic": [
            [0, 3],
            [0, 4],
            [0, 5],
            [1, 3],
            [1, 4],
            [1, 5],
            [2, 3],
            [2, 4],
            [2, 5],
            [3, 4],
            [3, 5],
            [4, 5],
        ],
        "hexagonal": [[0, 5], [1, 5], [2, 3], [2, 4], [2, 5], [3, 4]],
        "rhombohedral": [[0, 5], [1, 5], [2, 3], [2, 4], [2, 5], [3, 4]],
        "monoclinic": [[0, 4], [0, 5], [1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5]],
        "triclinic": [],
    }
    
    
    equal_elements: Dict[str, List[List[List[int]]]] = {
        "cubic": [
            [[0, 0], [1, 1]],
            [[0, 0], [2, 2]],
            [[1, 1], [2, 2]],
            [[0, 1], [0, 2]],
            [[0, 1], [1, 2]],
            [[0, 2], [1, 2]],
            [[3, 3], [4, 4]],
            [[3, 3], [5, 5]],
            [[4, 4], [5, 5]],
        ],
        "tetragonal": [
            [[0, 0], [1, 1]],
            [[0, 2], [1, 2]],
            [[3, 3], [4, 4]],
        ],
        "orthorhombic": [],
        "hexagonal": [
            [[0, 0], [1, 1]],
            [[0, 2], [1, 2]],
            [[3, 3], [4, 4]],
        ],
        "rhombohedral": [
            [[0, 0], [1, 1]],
            [[0, 2], [1, 2]],
            [[3, 3], [4, 4]],
            [[0, 3], [4, 5]],
            [[1, 4], [3, 5]],
        ],
        "monoclinic": [],
        "triclinic": [],
    }
    
    
    inverse_elements: Dict[str, List[List[List[int]]]] = {
        "cubic": [],
        "tetragonal": [[[0, 5], [1, 5]]],
        "orthorhombic": [],
        "hexagonal": [],
        "rhombohedral": [
            [[0, 3], [1, 3]],
            [[4, 5], [1, 3]],
            [[1, 4], [0, 4]],
            [[3, 5], [0, 4]],
        ],
        "monoclinic": [],
        "triclinic": [],
    }
    
