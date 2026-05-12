# matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients#

matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients(_time : ndarray_, _msd : ndarray_, _dimensions : int = 1_) → Tuple[float, float][[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#fit_diffusion_coefficients)#
    

Obtain the diffusion coefficient from fitting the time and the mean squared displacement.

Parameters
    

  * **time** (_np.ndarray_) – The time of each step in fs.

  * **msd** (_np.ndarray_) – The mean squared displacement.

  * **dimensions** (_int_ _,__optional_) – The dimensions. Defaults to 1.



Returns
    

The diffusion coefficient and its standard deviation.

Return type
    

Tuple[float, float]

[ __ previous matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult ](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html "previous page") [ next matlantis_features.features.md.post_features.diffusion.get_msd __](matlantis_features.features.md.post_features.diffusion.get_msd.html "next page")
