# matlantis_features.features.phonon.symmetry.get_permutations#

matlantis_features.features.phonon.symmetry.get_permutations(_pos_a : ndarray_, _pos_b : ndarray_, _symprec : float = 1e-05_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/symmetry.html#get_permutations)#
    

Get the permuation index which enables pos_b = pos_a[ind].

Parameters
    

  * **pos_a** (_np.ndarray_) – coordinations before permutation.

  * **pos_b** (_np.ndarray_) – coordinations after permutation.

  * **symprec** (_float_ _,__optional_) – distance tolerance to find crystal symmetry. Defaults to 1e-05.



Returns
    

the permuation index.

Return type
    

np.ndarray

[ __ previous matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry ](matlantis_features.features.phonon.symmetry.get_least_displacements_from_symmetry.html "previous page") [ next matlantis_features.features.phonon.symmetry.resolve_force_constants __](matlantis_features.features.phonon.symmetry.resolve_force_constants.html "next page")
