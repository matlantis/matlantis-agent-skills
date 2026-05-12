# matlantis_features.features.elasticity.tensors.StrainTensorList#

_class _matlantis_features.features.elasticity.tensors.StrainTensorList(_data : List[List[[StrainTensor](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor")]]_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensorList)#
    

Bases: `object`

Tensor list object for the strain tensors.

Methods

`__init__`(data) | Initialize an instance.  
---|---  
`from_strains`(diagonal_strains, ...) | Create StrainTensorList from the list of diagonal and off-diagonal strain values.  
`get_samples`(vix) | Get the strain sample list for the fitting calculation.  
  
__init__(_data : List[List[[StrainTensor](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor")]]_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensorList.__init__)#
    

Initialize an instance.

Parameters
    

**data** (_list_ _[__list_ _[_[_StrainTensor_](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor") _]__]_) – 2D-list of the strain tensors.

_classmethod _from_strains(_diagonal_strains : List[float]_, _off_diagonal_strains : List[float]_) → StrainTensorList[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensorList.from_strains)#
    

Create StrainTensorList from the list of diagonal and off-diagonal strain values.

Parameters
    

  * **diagonal_strains** (_list_ _[__float_ _]_) – list of diagonal strain values.

  * **off_diagonal_strains** (_list_ _[__float_ _]_) – list of diagonal strain values.



Returns
    

StrainTensorList object.

Return type
    

StrainTensorList

get_samples(_vix : int_) → List[float][[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensorList.get_samples)#
    

Get the strain sample list for the fitting calculation.

Parameters
    

**vix** (_int_) – Voigt index corresponding to the data-point axis.

Returns
    

List of stress values for the fitting calculation.

Return type
    

list[float]

[ __ previous matlantis_features.features.elasticity.tensors.StrainTensor ](matlantis_features.features.elasticity.tensors.StrainTensor.html "previous page") [ next matlantis_features.features.elasticity.tensors.StressTensor __](matlantis_features.features.elasticity.tensors.StressTensor.html "next page")
