# matlantis_features.features.common.vibration.VibrationFeatureResult#

_class _matlantis_features.features.common.vibration.VibrationFeatureResult(_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_, _indices : Optional[List[int]]_, _delta : float_, _nfree : int_, _method : str_, _direction : str_, _all_forces : ndarray_, _force_constant : ndarray_, _frequencies : ndarray_, _image_frequencies : ndarray_, _mode : ndarray_, _zero_point_energy : float_, _potential_energy : float_)[[source]](../_modules/matlantis_features/features/common/vibration.html#VibrationFeatureResult)#
    

Bases: `object`

A dataclass for result of VibrationFeature.

Methods

`__init__`(atoms, indices, delta, nfree, ...) |   
---|---  
  
__init__(_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_, _indices : Optional[List[int]]_, _delta : float_, _nfree : int_, _method : str_, _direction : str_, _all_forces : ndarray_, _force_constant : ndarray_, _frequencies : ndarray_, _image_frequencies : ndarray_, _mode : ndarray_, _zero_point_energy : float_, _potential_energy : float_) → None#
    

Attributes

dict#
    

Dictionary representation of the VibrationFeatureResult.

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
    

all_forces _: ndarray_#
    

force_constant _: ndarray_#
    

frequencies _: ndarray_#
    

image_frequencies _: ndarray_#
    

mode _: ndarray_#
    

zero_point_energy _: float_#
    

potential_energy _: float_#
    

[ __ previous matlantis_features.features.common.vibration.VibrationFeature ](matlantis_features.features.common.vibration.VibrationFeature.html "previous page") [ next matlantis_features.features.common.vibration.FastVibrations __](matlantis_features.features.common.vibration.FastVibrations.html "next page")
