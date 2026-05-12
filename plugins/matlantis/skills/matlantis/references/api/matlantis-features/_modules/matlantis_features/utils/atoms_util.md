# Source code for matlantis_features.utils.atoms_util
    
    
    from enum import Enum, unique
    
    import numpy as np
    from ase import Atoms
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.convert_atoms_to_upper.html#matlantis_features.utils.atoms_util.convert_atoms_to_upper)def convert_atoms_to_upper(atoms: Atoms) -> Atoms:
        """Rotate atoms
        ase NPT module requires `atoms` to have following upper-triangular form for `cell`:
        [[X X X]
        [0 X X]
        [0 0 X]]
        this method rotates `atoms` to make its cell upper triangular.
    
        Args:
            atoms (Atoms): input atoms.
        Returns:
            rotated_atoms (Atoms): output, rotated atoms.
        """
        if not np.all(atoms.get_pbc() == np.array([True, True, True])):
            raise ValueError("The method is only applicable for the structure with PBC.")
    
        rotated_atoms = atoms.copy()
    
        # cell "c" -> z-axis
        rotated_atoms.rotate(rotated_atoms.cell[2], (0, 0, 1), rotate_cell=True)
    
        # cell "b" -> yz-plane
        bx, by, bz = rotated_atoms.cell[1, :]
        angle = 90.0 - np.rad2deg(np.arctan2(by, bx))
        rotated_atoms.rotate(angle, "z", rotate_cell=True)
        # [Note] cell "a" can be arbitrary.
    
        # supress numerical precision, lower triangular values must be 0.0.
        # ase.atoms.rotate need values inside the vector to be rotated to be >= 1e-6
        rotated_atoms.cell = np.where(np.abs(rotated_atoms.cell) < 1e-6, 0.0, rotated_atoms.cell)
    
        # Check that cell is uppper triangular form.
        m = rotated_atoms.cell
        assert m[1, 0] == m[2, 0] == m[2, 1] == 0.0, f"cell {m} is not upper triangular!"
        return rotated_atoms
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.get_inertia_tensor.html#matlantis_features.utils.atoms_util.get_inertia_tensor)def get_inertia_tensor(atoms: Atoms) -> np.ndarray:
        """Calculate the inertia tensor of molecule.
    
        Args:
            atoms (Atoms): The input structure.
        Returns:
            np.ndarray : The inertia tensor.
        """
        com = atoms.get_center_of_mass()
        r = atoms.positions - com
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        m = atoms.get_masses()
        Ixx = np.sum(m * (y**2 + z**2))
        Iyy = np.sum(m * (x**2 + z**2))
        Izz = np.sum(m * (x**2 + y**2))
        Ixy = -np.sum(m * x * y)
        Iyz = -np.sum(m * y * z)
        Ixz = -np.sum(m * x * z)
        return np.array(  # type: ignore
            [
                [Ixx, Ixy, Ixz],
                [Ixy, Iyy, Iyz],
                [Ixz, Iyz, Izz],
            ]
        )
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.get_total_inertia.html#matlantis_features.utils.atoms_util.get_total_inertia)def get_total_inertia(atoms: Atoms) -> float:
        """Calculate the moment of inertia of molecule.
    
        Args:
            atoms (Atoms): The input structure.
        Returns:
            float : The momenta of inertia.
        """
        com = atoms.get_center_of_mass()
        r = np.linalg.norm(atoms.positions - com, axis=1)
        total_inertia: float = np.sum(atoms.get_masses() * r**2)
        return total_inertia
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.MoleculeLinearity.html#matlantis_features.utils.atoms_util.MoleculeLinearity)@unique
    class MoleculeLinearity(Enum):
        monatomic = 0
        linear = 1
        nonlinear = 2
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.check_linearity.html#matlantis_features.utils.atoms_util.check_linearity)def check_linearity(atoms: Atoms, eigen_tolerance: float = 0.0001) -> MoleculeLinearity:
        """Check the linearity of molecule. The linearity can be linear, nonlinear and monatomic.
    
        Args:
            atoms (Atoms): The input structure.
            eigen_tolerance (float, optional): Tolerance of eigen value of inertia tensor when
                counting the number of zero elements. Defaults to 0.0001.
        Returns:
            MoleculeLinearity : The linearity of the molecule.
        """
        if len(atoms) == 1:
            return MoleculeLinearity.monatomic
        inertia_tensor = get_inertia_tensor(atoms) / get_total_inertia(atoms)
        w, v = np.linalg.eigh(inertia_tensor)
        if np.sum(w < eigen_tolerance) == 1:
            return MoleculeLinearity.linear
        elif np.sum(w < eigen_tolerance) == 0:
            return MoleculeLinearity.nonlinear
        else:
            raise ValueError("Cannot identify the linearity of the molecule.")
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.atoms_util.check_symmetry_number.html#matlantis_features.utils.atoms_util.check_symmetry_number)def check_symmetry_number(atoms: Atoms) -> int:
        """Calculate the symmetry number of molecule.
    
        Args:
            atoms (Atoms): The input structure.
        Returns:
            int : The symmetry number.
        """
        _atoms = atoms.copy()  # copy to avoid overwriting user-defined magmoms
        _anums = _atoms.get_atomic_numbers()
        _magmom = np.zeros_like(_anums)
        if _anums.sum() % 2 == 1:
            _magmom[0] = 1
        _atoms.set_initial_magnetic_moments(_magmom)
    
        mol = AseAtomsAdaptor().get_molecule(_atoms)
        point_group = PointGroupAnalyzer(mol).get_pointgroup()
        symbol = point_group.__str__()
        if symbol in ["C1", "Ci", "Cs", "C*v"]:
            return 1
        elif symbol in ["D*h"]:
            return 2
        elif symbol[0] == "S":
            return int(symbol[1]) // 2
        elif symbol[0] == "C":
            return int(symbol[1])
        elif symbol[0] == "D":
            return int(symbol[1]) * 2
        elif symbol[0] == "T":
            return 12
        elif symbol in ["Oh"]:
            return 24
        elif symbol in ["Ih"]:
            return 60
        else:
            return 1
    
    
    
