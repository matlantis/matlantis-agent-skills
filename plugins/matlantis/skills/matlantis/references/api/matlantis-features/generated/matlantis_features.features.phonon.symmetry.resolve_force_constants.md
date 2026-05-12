# matlantis_features.features.phonon.symmetry.resolve_force_constants#

matlantis_features.features.phonon.symmetry.resolve_force_constants(_atoms : Atoms_, _atoms_N : Atoms_, _symmetry_operations : [SymmetryOperations](matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations "matlantis_features.features.phonon.symmetry.SymmetryOperations")_, _forces : ndarray_, _disp_index : ndarray_, _disp : ndarray_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#resolve_force_constants)#
    

Calculate the force constant from the least number of displacemented structures.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **atoms_N** (_Atoms_) – The supercell.

  * **symmetry_operations** ([_SymmetryOperations_](matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations "matlantis_features.features.phonon.symmetry.SymmetryOperations")) – The symmetry operations in the supercell.

  * **forces** (_np.ndarray_) – The forces of the displaced structures.

  * **disp_index** (_np.ndarray_) – The index of the displaced atom.

  * **disp** (_np.ndarray_) – the displacment vector of the displaced atom.



Returns
    

The force constant.

Return type
    

np.ndarray

[ __ previous matlantis_features.features.phonon.symmetry.get_permutations ](matlantis_features.features.phonon.symmetry.get_permutations.html "previous page") [ next matlantis_features.features.phonon.symmetry.similarity_transformation __](matlantis_features.features.phonon.symmetry.similarity_transformation.html "next page")
