# matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult#

_class _matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult(_volume : [VolumeResult](matlantis_features.features.md.post_features.thermal_expansion.VolumeResult.html#matlantis_features.features.md.post_features.thermal_expansion.VolumeResult "matlantis_features.features.md.post_features.thermal_expansion.VolumeResult")_, _thermal_expansion : [ThermalExpansionResult](matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult.html#matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult "matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult")_)[[source]](../_modules/matlantis_features/features/md/post_features/thermal_expansion.html#PostMDThermalExpansionFeatureResult)#
    

Bases: `object`

Result dataclass for the thermal expansion feature.

Methods

`__init__`(volume, thermal_expansion) |   
---|---  
`plot`([plt_name]) | Plot the result.  
  
__init__(_volume : [VolumeResult](matlantis_features.features.md.post_features.thermal_expansion.VolumeResult.html#matlantis_features.features.md.post_features.thermal_expansion.VolumeResult "matlantis_features.features.md.post_features.thermal_expansion.VolumeResult")_, _thermal_expansion : [ThermalExpansionResult](matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult.html#matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult "matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult")_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/thermal_expansion.html#PostMDThermalExpansionFeatureResult.plot)#
    

Plot the result.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

Attributes

volume _: [VolumeResult](matlantis_features.features.md.post_features.thermal_expansion.VolumeResult.html#matlantis_features.features.md.post_features.thermal_expansion.VolumeResult "matlantis_features.features.md.post_features.thermal_expansion.VolumeResult")_#
    

thermal_expansion _: [ThermalExpansionResult](matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult.html#matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult "matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult")_#
    

[ __ previous matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature ](matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature.html "previous page") [ next matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult __](matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult.html "next page")
