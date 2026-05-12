# matlantis_features.features.base.FeatureBaseCaller#

_class _matlantis_features.features.base.FeatureBaseCaller(_name : str_, _base : Any_, _attr : Any_)[[source]](../_modules/matlantis_features/features/base.html#FeatureBaseCaller)#
    

Bases: `type`

The metaclass used in constructing the feature.

Methods

`__init__`(*args, **kwargs) |   
---|---  
`__call__`(*args, **kwargs) | Call self as a function.  
`call_decorate`(derived_call) | Decorate the __call__ function of the feature.  
`mro`() | Return a type's method resolution order.  
  
__init__(_* args_, _** kwargs_)#
    

__call__(_* args_, _** kwargs_)#
    

Call self as a function.

_classmethod _call_decorate(_derived_call : Callable[[...], T]_) → Callable[[...], T][[source]](../_modules/matlantis_features/features/base.html#FeatureBaseCaller.call_decorate)#
    

Decorate the __call__ function of the feature.

Parameters
    

**derived_call** (_...__- > T_) – The original __call__ function of the feature.

Returns
    

The decorated __call__ function.

Return type
    

… -> T

mro()#
    

Return a type’s method resolution order.

[ __ previous matlantis_features.features.base.FeatureBase ](matlantis_features.features.base.FeatureBase.html "previous page") [ next matlantis_features.features.base.FeatureBaseOutput __](matlantis_features.features.base.FeatureBaseOutput.html "next page")
