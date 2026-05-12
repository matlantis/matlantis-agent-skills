# matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult#

_class _matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult(_diffusion_coefficient : Dict[str, Dict[str, float]]_, _diffusion_coefficient_std : Dict[str, Dict[str, float]]_, _mean_squared_displacement : Dict[str, ndarray]_, _timestep : float_, _init_time : float_, _stride : int_, _number_of_segments : int_)[[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#OldPostMDDiffusionFeatureResult)#
    

Bases: `object`

A dataclass for result of PostMDDiffusionFeature.

Methods

`__init__`(diffusion_coefficient, ...) |   
---|---  
`plot`([plt_name]) | Plot the result.  
  
__init__(_diffusion_coefficient : Dict[str, Dict[str, float]]_, _diffusion_coefficient_std : Dict[str, Dict[str, float]]_, _mean_squared_displacement : Dict[str, ndarray]_, _timestep : float_, _init_time : float_, _stride : int_, _number_of_segments : int_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#OldPostMDDiffusionFeatureResult.plot)#
    

Plot the result.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

Attributes

dict#
    

Convert the result to the dict object.

Returns
    

Resulting dict object.

Return type
    

dict[str, Any]

diffusion_coefficient _: Dict[str, Dict[str, float]]_#
    

diffusion_coefficient_std _: Dict[str, Dict[str, float]]_#
    

mean_squared_displacement _: Dict[str, ndarray]_#
    

timestep _: float_#
    

init_time _: float_#
    

stride _: int_#
    

number_of_segments _: int_#
    

[ __ previous matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature ](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature.html "previous page") [ next matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult __](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html "next page")
