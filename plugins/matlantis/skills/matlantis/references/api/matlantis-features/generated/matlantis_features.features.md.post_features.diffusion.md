# matlantis_features.features.md.post_features.diffusion#

PostFeatures

[`OldPostMDDiffusionFeature`](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature "matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature")() | The matlantis-feature for calculating the diffusion coefficient using the MD trajectory.  
---|---  
[`PostMDDiffusionFeature`](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature "matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature")() | The matlantis-feature for calculating the diffusion coefficient using the MD trajectory.  
  
FeatureResults

[`OldPostMDDiffusionFeatureResult`](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult "matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult")(...) | A dataclass for result of PostMDDiffusionFeature.  
---|---  
[`PostMDDiffusionFeatureResult`](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult "matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult")(...) | A dataclass for result of PostMDDiffusionFeature.  
  
Functions

[`fit_diffusion_coefficients`](matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients.html#matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients "matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients")(time, msd[, ...]) | Obtain the diffusion coefficient from fitting the time and the mean squared displacement.  
---|---  
[`get_msd`](matlantis_features.features.md.post_features.diffusion.get_msd.html#matlantis_features.features.md.post_features.diffusion.get_msd "matlantis_features.features.md.post_features.diffusion.get_msd")(positions_array, masses, element_indices) | Calculate the mean squared displacement.  
  
[ __ previous matlantis_features.features.md.md_system_base.MDSystemBase ](matlantis_features.features.md.md_system_base.MDSystemBase.html "previous page") [ next matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature __](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html "next page")
