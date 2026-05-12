# matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult#

_class _matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult(_num_of_segments : int_, _viscosity : Dict[str, List[float]]_, _viscosity_std : Dict[str, List[float]]_, _time : ndarray_, _acf_segment : ndarray_, _acf_avg : ndarray_, _vis_segment : ndarray_, _vis_avg : ndarray_)[[source]](../_modules/matlantis_features/features/md/post_features/emd_viscosity.html#PostEMDViscosityFeatureResult)#
    

Bases: `object`

A dataclass for result of PostEMDViscosityFeature.

Methods

`__init__`(num_of_segments, viscosity, ...) |   
---|---  
`plot`([plt_name]) | Plot the result.  
`to_dataframe`([csv_name]) | Convert the result to the Pandas DataFrame object.  
  
__init__(_num_of_segments : int_, _viscosity : Dict[str, List[float]]_, _viscosity_std : Dict[str, List[float]]_, _time : ndarray_, _acf_segment : ndarray_, _acf_avg : ndarray_, _vis_segment : ndarray_, _vis_avg : ndarray_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/emd_viscosity.html#PostEMDViscosityFeatureResult.plot)#
    

Plot the result.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/md/post_features/emd_viscosity.html#PostEMDViscosityFeatureResult.to_dataframe)#
    

Convert the result to the Pandas DataFrame object.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – CSV file name to write the DataFrame. If None, no CSV file is generated. Defaults to None.

Returns
    

Resulting Pandas DataFrame object.

Return type
    

pd.DataFrame

Attributes

num_of_segments _: int_#
    

viscosity _: Dict[str, List[float]]_#
    

viscosity_std _: Dict[str, List[float]]_#
    

time _: ndarray_#
    

acf_segment _: ndarray_#
    

acf_avg _: ndarray_#
    

vis_segment _: ndarray_#
    

vis_avg _: ndarray_#
    

[ __ previous matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature ](matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature.html "previous page") [ next matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity __](matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity.html "next page")
