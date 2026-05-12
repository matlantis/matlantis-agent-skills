# matlantis_features.features.base.FeatureBaseResult#

_class _matlantis_features.features.base.FeatureBaseResult(_feature : [FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")_, _call_params : Dict[str, Any]_, _output : [FeatureBaseOutput](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseResult)#
    

Bases: `object`

The base of the matlantis features result classes.

Methods

`__init__`(feature, call_params, output, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseResult object from serialized dict.  
`load`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`save`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseResult.  
  
__init__(_feature : [FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")_, _call_params : Dict[str, Any]_, _output : [FeatureBaseOutput](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_) → None#
    

_classmethod _from_dict(_res : Dict[str, Any]_) → FeatureBaseResult[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseResult.from_dict)#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBaseResult from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

FeatureBaseResult

_classmethod _load(_filename : str_) → FeatureBaseResult[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseResult.load)#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**filename** (_str_) – Name of file to save the FeatureBaseResult to.

Returns
    

The FeatureBaseResult object if loading was successful.

Return type
    

FeatureBaseResult

save(_filename : str_) → bool[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseResult.save)#
    

Construct a FeatureBaseResult object from serialized dict.

Parameters
    

**filename** (_str_) – Name of file to save the FeatureBaseResult to.

Returns
    

Whether saving the result was successful or not.

Return type
    

bool

to_dict() → Dict[str, Any][[source]](../_modules/matlantis_features/features/base.html#FeatureBaseResult.to_dict)#
    

Dictionary representation of the FeatureBaseResult.

Returns
    

A dict containing a serialized FeatureBaseResult.

Return type
    

dict[str, Any]

Attributes

feature _: [FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")_#
    

call_params _: Dict[str, Any]_#
    

output _: [FeatureBaseOutput](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")_#
    

calculator_info _: Optional[Dict[str, Any]]_#
    

version_info _: Dict[str, Any]_#
    

user_info _: Dict[str, Any]_#
    

[ __ previous matlantis_features.features.base.FeatureBaseOutput ](matlantis_features.features.base.FeatureBaseOutput.html "previous page") [ next matlantis_features.features.base.copy_to_child_conversion __](matlantis_features.features.base.copy_to_child_conversion.html "next page")
