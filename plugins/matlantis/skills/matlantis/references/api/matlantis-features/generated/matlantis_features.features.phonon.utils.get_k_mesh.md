# matlantis_features.features.phonon.utils.get_k_mesh#

matlantis_features.features.phonon.utils.get_k_mesh(_size : ndarray_, _method : str = 'mp'_, _shift : Optional[ndarray] = None_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/utils.html#get_k_mesh)#
    

Construct a uniform sampling of k-space of given size.

Parameters
    

  * **size** (_np.ndarray_) – The size of k-point mesh grid.

  * **method** (_str_ _,__optional_) – Create the k-point mesh grid with Monkhorst-Pack scheme (‘mp”) or Gamma centered scheme (‘gamma’). Defaults to ‘mp’.

  * **shift** (_np.ndarray_ _or_ _None_ _,__optional_) – An array of shape (3,) containing the shift of the mesh grid in the direction along the corresponding reciprocal axes. Defaults to None.



Returns
    

The k-point mesh.

Return type
    

np.ndarray

[ __ previous matlantis_features.features.phonon.utils.PhononFrequency ](matlantis_features.features.phonon.utils.PhononFrequency.html "previous page") [ next matlantis_features.features.phonon.utils.get_primitive_structure __](matlantis_features.features.phonon.utils.get_primitive_structure.html "next page")
