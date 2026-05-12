# matlantis_features.features.common.neb.NEBFeatureOutput#

_class _matlantis_features.features.common.neb.NEBFeatureOutput(_neb_images : List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _converged : bool_, _steps : int_, _energy_log : List[List[float]]_, _force_log : List[List[float]]_)[[source]](../_modules/matlantis_features/features/common/neb.html#NEBFeatureOutput)#
    

Bases: [`FeatureBaseOutput`](matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput "matlantis_features.features.base.FeatureBaseOutput")

A dataclass for result output of NEBFeature.

Methods

`__init__`(neb_images, converged, steps, ...) |   
---|---  
`from_dict`(res) | Construct a FeatureBaseOutput object from serialized dict.  
`to_dict`() | Dictionary representation of the FeatureBaseOutput.  
  
__init__(_neb_images : List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _converged : bool_, _steps : int_, _energy_log : List[List[float]]_, _force_log : List[List[float]]_) → None#
    

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

[ __ previous matlantis_features.features.common.neb.NEBFeatureResult ](matlantis_features.features.common.neb.NEBFeatureResult.html "previous page") [ next Gas Thermo Features __](../matlantis_features.features.common.gas_thermo.html "next page")
