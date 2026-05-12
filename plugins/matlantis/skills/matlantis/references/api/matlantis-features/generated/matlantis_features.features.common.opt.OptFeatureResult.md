# matlantis_features.features.common.opt.OptFeatureResult#

_class _matlantis_features.features.common.opt.OptFeatureResult(_feature : [OptFeatureBase](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")_, _call_params : Dict[str, Any]_, _output : [OptFeatureOutput](matlantis_features.features.common.opt.OptFeatureOutput.html#matlantis_features.features.common.opt.OptFeatureOutput "matlantis_features.features.common.opt.OptFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/common/opt.html#OptFeatureResult)#
    

Bases: [`FeatureBaseResult`](matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult "matlantis_features.features.base.FeatureBaseResult")

Class with the results for various optimizer features.

Methods

`__init__`(feature, call_params, output, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseResult object from serialized dict.  
`load`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`save`(filename) | Construct a FeatureBaseResult object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseResult.  
  
__init__(_feature : [OptFeatureBase](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")_, _call_params : Dict[str, Any]_, _output : [OptFeatureOutput](matlantis_features.features.common.opt.OptFeatureOutput.html#matlantis_features.features.common.opt.OptFeatureOutput "matlantis_features.features.common.opt.OptFeatureOutput")_, _calculator_info : Optional[Dict[str, Any]]_, _version_info : Dict[str, Any]_, _user_info : Dict[str, Any]_) → None#
    

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

Attributes

feature _: [OptFeatureBase](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")_#
    

output _: [OptFeatureOutput](matlantis_features.features.common.opt.OptFeatureOutput.html#matlantis_features.features.common.opt.OptFeatureOutput "matlantis_features.features.common.opt.OptFeatureOutput")_#
    

call_params _: Dict[str, Any]_#
    

calculator_info _: Optional[Dict[str, Any]]_#
    

version_info _: Dict[str, Any]_#
    

user_info _: Dict[str, Any]_#
    

[ __ previous matlantis_features.features.common.opt.MDMinASEOptFeature ](matlantis_features.features.common.opt.MDMinASEOptFeature.html "previous page") [ next matlantis_features.features.common.opt.OptFeatureBase __](matlantis_features.features.common.opt.OptFeatureBase.html "next page")
