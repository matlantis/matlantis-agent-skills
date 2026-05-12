# matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult#

_class _matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult(_viscosity : Dict[str, List[float]]_, _viscosity_std : Dict[str, List[float]]_, _slab_velocity_x : List[float]_, _slab_z : List[float]_)[[source]](../_modules/matlantis_features/features/md/post_features/nemd_viscosity.html#PostNEMDViscosityFeatureResult)#
    

Bases: `object`

A dataclass for result of PostNEMDViscosityFeature.

Methods

`__init__`(viscosity, viscosity_std, ...) |   
---|---  
`plot`([plt_name]) | Plot the result.  
`to_dataframe`([csv_name]) | Convert the result to the Pandas DataFrame object.  
  
__init__(_viscosity : Dict[str, List[float]]_, _viscosity_std : Dict[str, List[float]]_, _slab_velocity_x : List[float]_, _slab_z : List[float]_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/nemd_viscosity.html#PostNEMDViscosityFeatureResult.plot)#
    

Plot the result.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/md/post_features/nemd_viscosity.html#PostNEMDViscosityFeatureResult.to_dataframe)#
    

Convert the result to the Pandas DataFrame object.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – CSV file name to write the DataFrame. If None, no CSV file is generated. Defaults to None.

Returns
    

Resulting Pandas DataFrame object.

Return type
    

pd.DataFrame

Attributes

viscosity _: Dict[str, List[float]]_#
    

viscosity_std _: Dict[str, List[float]]_#
    

slab_velocity_x _: List[float]_#
    

slab_z _: List[float]_#
    

[ __ previous matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult ](matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult.html "previous page") [ next matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity __](matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity.html "next page")
