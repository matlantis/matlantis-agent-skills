# matlantis_features.features.common.neb.NEBFeatureResult#

_class _matlantis_features.features.common.neb.NEBFeatureResult(_feature : [NEBFeature](matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature "matlantis_features.features.common.neb.NEBFeature")_, _call_params : Dict[str, Any]_, _output : [NEBFeatureOutput](matlantis_features.features.common.neb.NEBFeatureOutput.html#matlantis_features.features.common.neb.NEBFeatureOutput "matlantis_features.features.common.neb.NEBFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureResult)#
    

Bases: [`FeatureBaseResult`](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")

A dataclass for result of NEBFeature.

Methods

`__init__`(feature, call_params, output, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseResult object from serialized dict.  
`load`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`plot`([plt_name]) | Plots the energy profile along the NEB pathway.  
`save`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`to_dataframe`([csv_name]) | Saves the NEB information to pandas.DataFrame.  
`to_dict`() | Dictionary representation of the FeatureBaseResult.  
`to_gif`([gif_name, rotation]) | Creates gif file to show the movement of atoms along the NEB pathway.  
`to_traj`(traj_name) | Creates ase.io.Trajectory to save the movement of atoms.  
  
__init__(_feature : [NEBFeature](matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature "matlantis_features.features.common.neb.NEBFeature")_, _call_params : Dict[str, Any]_, _output : [NEBFeatureOutput](matlantis_features.features.common.neb.NEBFeatureOutput.html#matlantis_features.features.common.neb.NEBFeatureOutput "matlantis_features.features.common.neb.NEBFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_) → None#
    

_classmethod _from_dict(_res : Dict[str, Any]_) → [FeatureBaseResult](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBaseResult from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

[FeatureBaseResult](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")

_classmethod _load(_filename : str_) → [FeatureBaseResult](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**filename** (_str_) – Name of file to save the FeatureBaseResult to.

Returns
    

The FeatureBaseResult object if loading was successful.

Return type
    

[FeatureBaseResult](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureResult.plot)#
    

Plots the energy profile along the NEB pathway.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of the energy profile will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

The plot of the NEB energy profile as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

save(_filename : str_) → bool#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**filename** (_str_) – Name of file to save the FeatureBaseResult to.

Returns
    

Whether saving the result was successful or not.

Return type
    

bool

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureResult.to_dataframe)#
    

Saves the NEB information to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the NEB information will be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains the energy and max force of each image during the optimization.

Return type
    

pd.DataFrame

to_dict() → Dict[str, Any]#
    

Dictionary representation of the FeatureBaseResult.

Returns
    

A dict containing a serialized FeatureBaseResult.

Return type
    

dict[str, Any]

to_gif(_gif_name : Optional[str] = None_, _rotation : str = '0x,0y,0z'_) → None[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureResult.to_gif)#
    

Creates gif file to show the movement of atoms along the NEB pathway.

Parameters
    

  * **gif_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the NEB pathway will be saved. Defaults to None.

  * **rotation** (_str_ _,__optional_) – The viewing angle. Defaults to ‘0x,0y,0z’.




to_traj(_traj_name : str_) → None[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureResult.to_traj)#
    

Creates ase.io.Trajectory to save the movement of atoms.

Parameters
    

**traj_name** (_str_) – The file name to which the NEB trajectory will be saved.

Attributes

feature _: [NEBFeature](matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature "matlantis_features.features.common.neb.NEBFeature")_#
    

output _: [NEBFeatureOutput](matlantis_features.features.common.neb.NEBFeatureOutput.html#matlantis_features.features.common.neb.NEBFeatureOutput "matlantis_features.features.common.neb.NEBFeatureOutput")_#
    

call_params _: Dict[str, Any]_#
    

calculator_info _: Optional[Dict[str, Any]]_#
    

version_info _: Dict[str, Any]_#
    

user_info _: Dict[str, Any]_#
    

[ __ previous matlantis_features.features.common.neb.NEBFeature ](matlantis_features.features.common.neb.NEBFeature.html "previous page") [ next matlantis_features.features.common.neb.NEBFeatureOutput __](matlantis_features.features.common.neb.NEBFeatureOutput.html "next page")
