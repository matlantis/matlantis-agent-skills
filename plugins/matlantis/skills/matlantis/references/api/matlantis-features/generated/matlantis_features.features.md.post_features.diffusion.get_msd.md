# matlantis_features.features.md.post_features.diffusion.get_msd#

matlantis_features.features.md.post_features.diffusion.get_msd(_positions_array : ndarray_, _masses : ndarray_, _element_indices : OrderedDict_, _com_indices : Optional[List[List[int]]] = None_, _direction : Optional[ndarray] = None_, _method : str = 'normal'_) → Tuple[ndarray, Optional[ndarray]][[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#get_msd)#
    

Calculate the mean squared displacement.

Parameters
    

  * **positions_array** (_np.ndarray_) – The atomic position of each atom in each step.

  * **masses** (_np.ndarray_) – The atomic number of each atom.

  * **element_indices** (_OrderedDict_) – The indices of atoms of each element.

  * **com_indices** (_Optional_ _[__List_ _[__List_ _[__int_ _]__]__]__,__optional_) – The indices of atoms who form a molecule. If it is used, the diffusion coefficient of the center of mass of molecules will be calculated. Defaults to None.

  * **direction** (_Optional_ _[__np.ndarray_ _]__,__optional_) – Calculate the diffusion coefficient along specific directions. It must be a nx3 array where n is number of directions. If None, only the diffusion coefficient along x, y and z directions is calculated. Defaults to None

  * **method** (_str_ _,__optional_) – The method to calculate mean square displacement (MSD). Now, two methods, i.e. “normal” and “segment”, are supported. Defaults to “normal”.



Returns
    

the mean square displacement of each atom, the mean square displacement of each molecule.

Return type
    

tuple[np.ndarray, np.ndarray or None]

[ __ previous matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients ](matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients.html "previous page") [ next matlantis_features.features.md.post_features.emd_viscosity __](matlantis_features.features.md.post_features.emd_viscosity.html "next page")
