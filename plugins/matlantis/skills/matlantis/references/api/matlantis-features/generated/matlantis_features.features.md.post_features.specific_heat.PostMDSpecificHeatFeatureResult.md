# matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult#

_class _matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult(_enthalpy : [EnthalpyResult](matlantis_features.features.md.post_features.specific_heat.EnthalpyResult.html#matlantis_features.features.md.post_features.specific_heat.EnthalpyResult "matlantis_features.features.md.post_features.specific_heat.EnthalpyResult")_, _specific_heat : [SpecificHeatResult](matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult.html#matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult "matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult")_)[[source]](../_modules/matlantis_features/features/md/post_features/specific_heat.html#PostMDSpecificHeatFeatureResult)#
    

Bases: `object`

Result dataclass for the specific heat and enthalpy.

Methods

`__init__`(enthalpy, specific_heat) |   
---|---  
`plot`([plt_name]) | Plot the result.  
  
__init__(_enthalpy : [EnthalpyResult](matlantis_features.features.md.post_features.specific_heat.EnthalpyResult.html#matlantis_features.features.md.post_features.specific_heat.EnthalpyResult "matlantis_features.features.md.post_features.specific_heat.EnthalpyResult")_, _specific_heat : [SpecificHeatResult](matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult.html#matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult "matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult")_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/specific_heat.html#PostMDSpecificHeatFeatureResult.plot)#
    

Plot the result.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

Attributes

enthalpy _: [EnthalpyResult](matlantis_features.features.md.post_features.specific_heat.EnthalpyResult.html#matlantis_features.features.md.post_features.specific_heat.EnthalpyResult "matlantis_features.features.md.post_features.specific_heat.EnthalpyResult")_#
    

specific_heat _: [SpecificHeatResult](matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult.html#matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult "matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult")_#
    

[ __ previous matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeature ](matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeature.html "previous page") [ next matlantis_features.features.md.post_features.specific_heat.EnthalpyResult __](matlantis_features.features.md.post_features.specific_heat.EnthalpyResult.html "next page")
