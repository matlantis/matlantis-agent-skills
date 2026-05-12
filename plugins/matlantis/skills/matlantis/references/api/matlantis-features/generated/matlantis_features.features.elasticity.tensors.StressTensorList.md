# matlantis_features.features.elasticity.tensors.StressTensorList#

_class _matlantis_features.features.elasticity.tensors.StressTensorList(_data : List[List[[StressTensor](matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor "matlantis_features.features.elasticity.tensors.StressTensor")]]_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StressTensorList)#
    

Bases: `object`

Tensor list object for the stress.

Methods

`__init__`(data) | Initialize an instance.  
---|---  
`get_samples`(vi, vj) | Get the stress sample list for the fitting calculation.  
  
__init__(_data : List[List[[StressTensor](matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor "matlantis_features.features.elasticity.tensors.StressTensor")]]_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StressTensorList.__init__)#
    

Initialize an instance.

Parameters
    

**data** (_list_ _[__list_ _[_[_StressTensor_](matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor "matlantis_features.features.elasticity.tensors.StressTensor") _]__]_) – 2D-list of the stress tensors.

get_samples(_vi : int_, _vj : int_) → List[float][[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StressTensorList.get_samples)#
    

Get the stress sample list for the fitting calculation.

Parameters
    

  * **vi** (_int_) – Voigt index corresponding to the data-point axis.

  * **vj** (_int_) – Voigt index corresponding to the stress tensor axis.



Returns
    

List of stress values for the fitting calculation.

Return type
    

list[float]

[ __ previous matlantis_features.features.elasticity.tensors.StressTensor ](matlantis_features.features.elasticity.tensors.StressTensor.html "previous page") [ next matlantis_features.features.elasticity.utils __](matlantis_features.features.elasticity.utils.html "next page")
