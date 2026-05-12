# matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput#

_class _matlantis_features.features.reaction.reaction_string.ReactionStringFeatureOutput(_num_segments : int_, _reaction_string_images : List[List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]]_, _transition_state_indices : List[Optional[int]]_, _imaginary_eigenvectors : List[ndarray]_, _converged : bool_)[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeatureOutput)#
    

Bases: [`FeatureBaseOutput`](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")

A dataclass containing the output of ReactionStringFeature.

Methods

`__init__`(num_segments, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseOutput object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseOutput.  
  
__init__(_num_segments : int_, _reaction_string_images : List[List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]]_, _transition_state_indices : List[Optional[int]]_, _imaginary_eigenvectors : List[ndarray]_, _converged : bool_) → None#
    

_classmethod _from_dict(_res : Dict[str, Any]_) → [FeatureBaseOutput](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")#
    

Construct a FeatureBaseOutput object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBaseOutput from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

[FeatureBaseOutput](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")

to_dict() → Dict[str, Any]#
    

Dictionary representation of the FeatureBaseOutput.

Returns
    

A dict containing a serialized FeatureBaseOutput.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult ](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html "previous page") [ next RestScan Features __](../matlantis_features.features.reaction.rest_scan.html "next page")
