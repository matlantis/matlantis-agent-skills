# matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix#

matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix(_atoms : Atoms_, _primitive_matrix : ndarray_, _prec : float = 0.001_) → Tuple[Atoms, ndarray][[source]](../_modules/matlantis_features/features/phonon/utils.html#get_primitive_structure_from_primitive_matrix)#
    

Transforms the input unit cell to the primitive cell according to the axis provided.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **primitive_matrix** (_np.ndarray_) – the primitive cell basis vectors.

  * **prec** (_float_ _,__optional_) – Wrapper the input structure into the primitive cell with certain precision. Defaults to 0.001.



Returns
    

The primitive cell.

Return type
    

tuple[Atoms, np.ndarray]

[ __ previous matlantis_features.features.phonon.utils.get_primitive_structure_auto ](matlantis_features.features.phonon.utils.get_primitive_structure_auto.html "previous page") [ next NEB Features __](../matlantis_features.features.common.neb.html "next page")
