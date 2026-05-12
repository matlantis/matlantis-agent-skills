# matlantis_features.atoms#

_class _matlantis_features.atoms.MatlantisAtoms(_atoms : Atoms_)[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms)#
    

Bases: `object`

MatlantisAtoms object is a collection of atoms used in matlantis-features and wraps ase.Atoms class.

__init__(_atoms : Atoms_)[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.__init__)#
    

Initialize an instance.

Parameters
    

**atoms** (_ASEAtoms_) – A collection of atoms specified by AseAtoms.

_property _ase_atoms _: Atoms_#
    

Get the AseAtoms corresponding to the object.

Returns
    

The AseAtoms corresponding to the object.

Return type
    

ASEAtoms

copy() → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.copy)#
    

Get a copy of the MatlantisAtoms object.

Returns
    

Copy of the atoms.

Return type
    

MatlantisAtoms

_classmethod _from_ase_atoms(_atoms : Atoms_) → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.from_ase_atoms)#
    

A classmethod to convert ASEAtoms to MatlantisAtoms.

Parameters
    

**atoms** (_ASEAtoms_) – AseAtoms object to convert.

Returns
    

MatlantisAtoms corresponding to the input AseAtoms.

Return type
    

MatlantisAtoms

_classmethod _from_dict(_res : Dict[str, Any]_) → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.from_dict)#
    

Construct a MatlantisAtoms object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized MatlantisAtoms from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

MatlantisAtoms

_classmethod _from_file(_file : Union[str, Path]_) → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.from_file)#
    

A classmethod to get MatlantisAtoms from a file using ase.io.read function.

Parameters
    

**file** (_str_ _or_ _pathlib.Path_) – A path to file to convert.

Returns
    

MatlantisAtoms corresponding to the input file.

Return type
    

MatlantisAtoms

_classmethod _from_rdkit_mol(_mol : Mol_) → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.from_rdkit_mol)#
    

A classmethod to get MatlantisAtoms from an rdkit.Chem.Mol object. The positions of atoms should be set.

Parameters
    

**mol** (_Mol_) – An rdkit.Chem.Mol object to convert.

Returns
    

MatlantisAtoms corresponding to the input Mol object.

Return type
    

MatlantisAtoms

_classmethod _from_smiles(_smiles : str_, _method : str = 'etkdg_v2'_, _seed : int = -1_) → MatlantisAtoms[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.from_smiles)#
    

A classmethod to get MatlantisAtoms from a SMILES.

Parameters
    

  * **smiles** (_str_) – A SMILES string to convert.

  * **method** (_str_ _,__optional_) – An algorithm for position embedding. The algorithm can be selected from `["etkdg_v2", "etkdg", "etdg", "kdg", "2d"]`. Defaults to ‘etkdg_v2’.

  * **seed** (_int_ _,__optional_) – A random seed for position embedding. Defaults to (-1).



Returns
    

MatlantisAtoms corresponding to the input SMILES.

Return type
    

MatlantisAtoms

_property _numbers _: ndarray_#
    

Get the array of atomic numbers of the atoms.

Returns
    

The array of atomic numbers of the atoms.

Return type
    

np.ndarray

_property _positions _: ndarray_#
    

Get the array of coordinates of the atoms.

Returns
    

The array of coordinates of the atoms.

Return type
    

np.ndarray

rotate_atoms_to_upper() → None[[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.rotate_atoms_to_upper)#
    

Rotate the MatlantisAtoms inplace to make its cell a upper triangular matrix.

to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/atoms.html#MatlantisAtoms.to_dict)#
    

Dictionary representation of the MatlantisAtoms.

Returns
    

The dictionary containing the serialized MatlantisAtoms.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.reaction.rest_scan.RestScanFeatureResult ](generated/matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html "previous page") [ next matlantis_features.filters __](matlantis_features.filters.html "next page")
