# matlantis_features.features.phonon.utils.get_primitive_structure#

matlantis_features.features.phonon.utils.get_primitive_structure(_atoms : Atoms_, _primitive_matrix : Union[ndarray, str]_, _prec : float = 0.001_) → Atoms[[source]](../_modules/matlantis_features/features/phonon/utils.html#get_primitive_structure)#
    

Transforms the input unit cell to the primitive cell.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **primitive_matrix** (_Union_ _[__np.ndarray_ _,__str_ _]_) – the primitive cell basis vectors or ‘auto’. If a 3x3 array is provided, the structure will be wrapped into the primitive cell according to the user-defined axis. If ‘auto’, the primitive cell will be automatically defined according to its symmetry.

  * **prec** (_float_ _,__optional_) – Wrapper the input structure into the primitive cell with certain precision. Defaults to 0.001.



Returns
    

The primitive cell.

Return type
    

Atoms

[ __ previous matlantis_features.features.phonon.utils.get_k_mesh ](matlantis_features.features.phonon.utils.get_k_mesh.html "previous page") [ next matlantis_features.features.phonon.utils.get_primitive_structure_auto __](matlantis_features.features.phonon.utils.get_primitive_structure_auto.html "next page")
