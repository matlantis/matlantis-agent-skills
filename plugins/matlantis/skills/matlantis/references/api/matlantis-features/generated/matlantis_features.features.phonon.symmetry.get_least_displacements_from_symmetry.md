# matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry#

matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry(_atoms : Atoms_, _atoms_N : Atoms_, _symmetry_operations : [SymmetryOperations](matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations "matlantis_features.features.phonon.symmetry.SymmetryOperations")_, _delta : float_) → Tuple[List[Atoms], ndarray, ndarray][[source]](../_modules/matlantis_features/features/phonon/symmetry.html#get_least_displacements_from_symmetry)#
    

To find the least number of displacements that need to calculate force constant.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **atoms_N** (_Atoms_) – The supercell.

  * **symmetry_operations** ([_SymmetryOperations_](matlantis_features.features.phonon.symmetry.SymmetryOperations.html#matlantis_features.features.phonon.symmetry.SymmetryOperations "matlantis_features.features.phonon.symmetry.SymmetryOperations")) – The symmetry operations in the supercell.

  * **delta** (_float_) – The distance of finite displacement applied to each atom.



Returns
    

A list of displaced structures, the index of the displaced atom, and the displacment vector of the displaced atom.

Return type
    

tuple[list[Atoms], np.ndarray, np.ndarray]

[ __ previous matlantis_features.features.phonon.symmetry.SymmetryOperations ](matlantis_features.features.phonon.symmetry.SymmetryOperations.html "previous page") [ next matlantis_features.features.phonon.symmetry.get_permutations __](matlantis_features.features.phonon.symmetry.get_permutations.html "next page")
