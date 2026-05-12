# matlantis_features.features.base.FeatureBaseOutput#

_class _matlantis_features.features.base.FeatureBaseOutput[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseOutput)#
    

Bases: `object`

The base of the matlantis features result output classes.

Methods

`__init__`() |   
---|---  
`from_dict`(res) | Construct a FeatureBaseOutput object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseOutput.  
  
__init__() → None#
    

_classmethod _from_dict(_res : Dict[str, Any]_) → FeatureBaseOutput[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseOutput.from_dict)#
    

Construct a FeatureBaseOutput object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBaseOutput from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

FeatureBaseOutput

to_dict() → Dict[str, Any][[source]](../_modules/matlantis_features/features/base.html#FeatureBaseOutput.to_dict)#
    

Dictionary representation of the FeatureBaseOutput.

Returns
    

A dict containing a serialized FeatureBaseOutput.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.base.FeatureBaseCaller ](matlantis_features.features.base.FeatureBaseCaller.html "previous page") [ next matlantis_features.features.base.FeatureBaseResult __](matlantis_features.features.base.FeatureBaseResult.html "next page")
