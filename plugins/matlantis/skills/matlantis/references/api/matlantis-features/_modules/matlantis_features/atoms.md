# Source code for matlantis_features.atoms
    
    
    import pathlib
    from dataclasses import dataclass
    from importlib import import_module
    from typing import Any, Dict, Union
    
    import ase
    import numpy as np
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.constraints import (
        ExternalForce,
        FixAtoms,
        FixBondLengths,
        FixCom,
        FixedLine,
        FixedMode,
        FixedPlane,
        FixInternals,
        FixLinearTriatomic,
        Hookean,
    )
    from ase.io import read
    from rdkit import Chem
    from rdkit.Chem import AllChem, Mol
    
    from matlantis_features.utils.atoms_util import convert_atoms_to_upper
    
    # TODO: fix
    ASEAtoms = Union[ase.Atoms]
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms)@dataclass
    class MatlantisAtoms(object):
        """MatlantisAtoms object is a collection of atoms used in matlantis-features \
        and wraps ase.Atoms class."""
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.__init__)    def __init__(self, atoms: ASEAtoms):
            """Initialize an instance.
    
            Args:
                atoms (ASEAtoms): A collection of atoms specified by AseAtoms.
            """
            self._atoms = atoms
    
    
    
        def __len__(self) -> int:
            return len(self._atoms)
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.from_ase_atoms)    @classmethod
        def from_ase_atoms(cls, atoms: ASEAtoms) -> "MatlantisAtoms":
            """A classmethod to convert ASEAtoms to MatlantisAtoms.
    
            Args:
                atoms (ASEAtoms): AseAtoms object to convert.
            Returns:
                MatlantisAtoms : MatlantisAtoms corresponding to the input AseAtoms.
            """
            return cls(atoms=atoms)
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.from_smiles)    @classmethod
        def from_smiles(cls, smiles: str, method: str = "etkdg_v2", seed: int = -1) -> "MatlantisAtoms":
            """A classmethod to get MatlantisAtoms from a SMILES.
    
            Args:
                smiles (str): A SMILES string to convert.
                method (str, optional):
                    An algorithm for position embedding. The algorithm can be selected from\
                    ``["etkdg_v2", "etkdg", "etdg", "kdg", "2d"]``. Defaults to 'etkdg_v2'.
                seed (int, optional): A random seed for position embedding. Defaults to (-1).
            Returns:
                MatlantisAtoms : MatlantisAtoms corresponding to the input SMILES.
            """
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
    
            if method not in ["etkdg_v2", "etkdg", "etdg", "kdg", "2d"]:
                raise ValueError("[ERROR] Unexpected value method={}".format(method))
    
            if method == "2d":
                AllChem.Compute2DCoords(mol)
            else:
                if method == "kdg":
                    params = AllChem.KDG()
                elif method == "etdg":
                    params = AllChem.ETDG()
                elif method == "etkdg":
                    params = AllChem.ETKDG()
                elif method == "etkdg_v2":
                    params = AllChem.ETKDGv2()
                params.randomSeed = seed
                AllChem.EmbedMolecule(mol, params)
    
            return cls.from_rdkit_mol(mol)
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.from_rdkit_mol)    @classmethod
        def from_rdkit_mol(cls, mol: Mol) -> "MatlantisAtoms":
            """A classmethod to get MatlantisAtoms from an rdkit.Chem.Mol object.\
            The positions of atoms should be set.
    
            Args:
                mol (Mol): An rdkit.Chem.Mol object to convert.
            Returns:
                MatlantisAtoms : MatlantisAtoms corresponding to the input Mol object.
            """
            symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
            positions = mol.GetConformer().GetPositions()
            atoms = Atoms(symbols, positions=positions)
            return cls(atoms=atoms)
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.from_file)    @classmethod
        def from_file(cls, file: Union[str, pathlib.Path]) -> "MatlantisAtoms":
            """A classmethod to get MatlantisAtoms from a file using ase.io.read function.
    
            Args:
                file (str or pathlib.Path): A path to file to convert.
            Returns:
                MatlantisAtoms : MatlantisAtoms corresponding to the input file.
            """
            atoms = read(file)
            return cls(atoms=atoms)
    
    
    
        @property
        def ase_atoms(self) -> ASEAtoms:
            """Get the AseAtoms corresponding to the object.
    
            Returns:
                ASEAtoms : The AseAtoms corresponding to the object.
            """
            return self._atoms
    
        @property
        def positions(self) -> np.ndarray:
            """Get the array of coordinates of the atoms.
    
            Returns:
                np.ndarray : The array of coordinates of the atoms.
            """
            return self._atoms.positions  # type:ignore
    
        @property
        def numbers(self) -> np.ndarray:
            """Get the array of atomic numbers of the atoms.
    
            Returns:
                np.ndarray : The array of atomic numbers of the atoms.
            """
            return self._atoms.numbers  # type:ignore
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.rotate_atoms_to_upper)    def rotate_atoms_to_upper(self) -> None:
            """Rotate the MatlantisAtoms inplace to make its cell a upper triangular matrix."""
            self._atoms = convert_atoms_to_upper(self.ase_atoms)
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.copy)    def copy(self) -> "MatlantisAtoms":
            """Get a copy of the MatlantisAtoms object.
    
            Returns:
                MatlantisAtoms: Copy of the atoms.
            """
            ase_atoms = self.ase_atoms.copy()
            ase_atoms.calc = None
            return self.__class__(ase_atoms)
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the MatlantisAtoms.
    
            Returns:
                dict[str, Any] : The dictionary containing the serialized MatlantisAtoms.
            """
            res = {}
            res = self.ase_atoms.copy().todict()
            res["cls_path"] = f"{self.__class__.__module__}.{self.__class__.__name__}"
            if hasattr(self.ase_atoms, "calc") and hasattr(self.ase_atoms.calc, "results"):
                res["results"] = self.ase_atoms.calc.results
            if "constraints" in res:
                for idx, constraint in enumerate(res["constraints"]):
                    if isinstance(
                        constraint,
                        (
                            ExternalForce,
                            FixAtoms,
                            FixBondLengths,
                            FixCom,
                            FixedLine,
                            FixedMode,
                            FixedPlane,
                            FixInternals,
                            FixLinearTriatomic,
                            Hookean,
                        ),
                    ):
                        res["constraints"][idx] = {
                            "cls_path": f"{constraint.__class__.__module__}.{constraint.__class__.__name__}",
                            **constraint.todict()["kwargs"],
                        }
                    elif hasattr(constraint, "to_dict"):  # For matlantis_features.symmetry.FixSymmetry, etc.
                        res["constraints"][idx] = constraint.to_dict()
            return res
    
    
    
    
    
    [[docs]](../../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms.from_dict)    @classmethod
        def from_dict(cls, res: Dict[str, Any]) -> "MatlantisAtoms":
            """Construct a MatlantisAtoms object from serialized dict.
    
            Args:
                res (dict[str, Any]): A dict containing a serialized MatlantisAtoms from to_dict().
            Returns:
                MatlantisAtoms : The deserialized object from provided dict.
            """
            ase_atoms = Atoms(positions=res["positions"], cell=res["cell"], numbers=res["numbers"], pbc=res["pbc"])
            if "results" in res:
                calc = SinglePointCalculator(ase_atoms)
                ase_atoms.calc = calc
                ase_atoms.calc.results = res["results"]
            if "constraints" in res:
                deserialized_constraints = []
                for constraint in res["constraints"].copy():
                    # if loaded by result.load(), constraints are no longer dict
                    # but if loaded directly from atoms.to_dict(), constraints are dict
                    if isinstance(constraint, dict):
                        module_name, cls_name = constraint["cls_path"].rsplit(".", 1)
                        if cls_name == "FixSymmetry":
                            constraint["atoms"] = cls.from_dict(constraint["atoms"])
                        del constraint["cls_path"]
                        Constraint = getattr(import_module(module_name), cls_name)
                        deserialized_constraints.append(Constraint(**constraint))
                    else:
                        deserialized_constraints.append(constraint)
                ase_atoms.constraints = deserialized_constraints
            return cls(ase_atoms)
    
    
    
