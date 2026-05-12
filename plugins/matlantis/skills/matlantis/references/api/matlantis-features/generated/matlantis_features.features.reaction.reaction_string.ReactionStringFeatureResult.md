# matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult#

_class _matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult(_feature : [ReactionStringFeature](matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature "matlantis_features.features.reaction.reaction_string.ReactionStringFeature")_, _call_params : Dict[str, Any]_, _output : [ReactionStringFeatureOutput](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput "matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeatureResult)#
    

Bases: [`FeatureBaseResult`](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")

A dataclass containing the result for ReactionStringFeature.

Methods

`__init__`(feature, call_params, output, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseResult object from serialized dict.  
`load`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`plot_energy`([segment_id, plt_name]) | Plot the energy profile along the reaction path.  
`save`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseResult.  
`to_gif`([gif_name, rotation]) | Creates gif file to show the movement of atoms along the reaction path.  
`to_traj`(traj_name) | Creates ase.io.Trajectory to save the reaction path.  
  
__init__(_feature : [ReactionStringFeature](matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature "matlantis_features.features.reaction.reaction_string.ReactionStringFeature")_, _call_params : Dict[str, Any]_, _output : [ReactionStringFeatureOutput](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput "matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_) → None#
    

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

plot_energy(_segment_id : Optional[int] = None_, _plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeatureResult.plot_energy)#
    

Plot the energy profile along the reaction path.

Parameters
    

  * **segment_id** (_int_ _or_ _None_ _,__optional_) – The index of the segment of path to plot. If None, all the path segments will be plotted togather. Defaults to None.

  * **plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of the energy profile will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.



Returns
    

The plot of the energy profile in the reaction string calculation result as a plotly.graph_objects.Figure instance.

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

to_dict() → Dict[str, Any]#
    

Dictionary representation of the FeatureBaseResult.

Returns
    

A dict containing a serialized FeatureBaseResult.

Return type
    

dict[str, Any]

to_gif(_gif_name : Optional[str] = None_, _rotation : str = '0x,0y,0z'_) → None[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeatureResult.to_gif)#
    

Creates gif file to show the movement of atoms along the reaction path.

Parameters
    

  * **gif_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the reaction path will be saved. Defaults to None.

  * **rotation** (_str_ _,__optional_) – The viewing angle. Defaults to ‘0x,0y,0z’.




to_traj(_traj_name : str_) → None[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeatureResult.to_traj)#
    

Creates ase.io.Trajectory to save the reaction path.

Parameters
    

**traj_name** (_str_) – The file name to which the NEB trajectory will be saved.

Attributes

feature _: [ReactionStringFeature](matlantis_features.features.reaction.reaction_string.ReactionStringFeature.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeature "matlantis_features.features.reaction.reaction_string.ReactionStringFeature")_#
    

output _: [ReactionStringFeatureOutput](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput "matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput")_#
    

call_params _: Dict[str, Any]_#
    

calculator_info _: Optional[Dict[str, Any]]_#
    

version_info _: Dict[str, Any]_#
    

user_info _: Dict[str, Any]_#
    

[ __ previous matlantis_features.features.reaction.reaction_string ](matlantis_features.features.reaction.reaction_string.html "previous page") [ next matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput __](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput.html "next page")
