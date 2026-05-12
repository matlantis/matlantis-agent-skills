# matlantis_features.features.common.hessian.HessianFeatureResult#

_class _matlantis_features.features.common.hessian.HessianFeatureResult(_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_, _indices : Optional[List[int]]_, _delta : float_, _nfree : int_, _method : str_, _direction : str_, _H : ndarray_)[[source]](../_modules/matlantis_features/features/common/hessian.html#HessianFeatureResult)#
    

Bases: `object`

A dataclass for result of HessianFeature.

Methods

`__init__`(atoms, indices, delta, nfree, ...) |   
---|---  
  
__init__(_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_, _indices : Optional[List[int]]_, _delta : float_, _nfree : int_, _method : str_, _direction : str_, _H : ndarray_) → None#
    

Attributes

dict#
    

Dictionary representation of the HessianFeatureResult.

Returns
    

The vibration calculation results in dictionary format.

Return type
    

dict[str, Any]

atoms _: [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_#
    

indices _: Optional[List[int]]_#
    

delta _: float_#
    

nfree _: int_#
    

method _: str_#
    

direction _: str_#
    

H _: ndarray_#
    

[ __ previous matlantis_features.features.common.hessian.HessianFeature ](matlantis_features.features.common.hessian.HessianFeature.html "previous page") [ next Phonon Features __](../matlantis_features.features.phonon.html "next page")
