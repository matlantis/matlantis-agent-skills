# matlantis_features.features.phonon.symmetry.SymmetryOperations#

_class _matlantis_features.features.phonon.symmetry.SymmetryOperations(_atoms : Atoms_, _symprec : float = 1e-05_)[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations)#
    

Bases: `object`

All the symmetry operations in the input structure.

Methods

`__init__`(atoms[, symprec]) | Initialize an instance.  
---|---  
`get_permutations`(index) | If the initial structure R is transformed with rotation[index] and translation[index] and get the rotated structure R'.  
`get_permutations_on_atom`(index) | If the initial strucutre R = {R_i} was rotated with matrix {r} on the atoms[index] and get the rotated structure R' = {R_i}', we can get a serial of index from this function that makes R'[ind] = R.  
`get_rotation_on_atom`(index[, return_indices]) | Find the rotation matrix center on the given atom.  
`get_rotations_coord`(index) | The rotation matrix for the cartesian coordinations instead of fractional coordinations.  
`get_translations_on_atom`(index) | To find a set of symmetry operation (both rotation and translation) that does not change the position of atoms[index].  
  
__init__(_atoms : Atoms_, _symprec : float = 1e-05_) → None[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.__init__)#
    

Initialize an instance.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **symprec** (_float_ _,__optional_) – Distance tolerance to find crystal symmetry. Defaults to 1e-05.




get_permutations(_index : int_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.get_permutations)#
    

If the initial structure R is transformed with rotation[index] and translation[index] and get the rotated structure R’. we can get a serial of index from this function which enable R’ = R[ind].

Parameters
    

**index** (_int_) – the index of symmetry operation.

Returns
    

the permutations index.

Return type
    

np.ndarray

get_permutations_on_atom(_index : int_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.get_permutations_on_atom)#
    

If the initial strucutre R = {R_i} was rotated with matrix {r} on the atoms[index] and get the rotated structure R’ = {R_i}’, we can get a serial of index from this function that makes R’[ind] = R.

Parameters
    

**index** (_int_) – the index of the center atom.

Returns
    

the permutations index.

Return type
    

np.ndarray

get_rotation_on_atom(_index : int_, _return_indices : bool = False_) → Union[ndarray, Tuple[ndarray, ndarray]][[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.get_rotation_on_atom)#
    

Find the rotation matrix center on the given atom.

Parameters
    

  * **index** (_int_) – The index of the center atom.

  * **return_indices** (_bool_ _,__optional_) – If True, also return the indices of rotation matrix. Defaults to False.



Returns
    

The rotation matrix center on the given atom.

Return type
    

np.ndarray or tuple[np.ndarray, np.ndarray]

get_rotations_coord(_index : int_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.get_rotations_coord)#
    

The rotation matrix for the cartesian coordinations instead of fractional coordinations.

Parameters
    

**index** (_int_) – The index of rotation operation.

Returns
    

The rotation matrix for the cartesian coordinations.

Return type
    

np.ndarray

get_translations_on_atom(_index : int_) → int[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#SymmetryOperations.get_translations_on_atom)#
    

To find a set of symmetry operation (both rotation and translation) that does not change the position of atoms[index].

Parameters
    

**index** (_int_) – the index of the unmoved atom

Returns
    

the index of rotation and translation operations.

Return type
    

int

[ __ previous matlantis_features.features.phonon.symmetry ](matlantis_features.features.phonon.symmetry.html "previous page") [ next matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry __](matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry.html "next page")
