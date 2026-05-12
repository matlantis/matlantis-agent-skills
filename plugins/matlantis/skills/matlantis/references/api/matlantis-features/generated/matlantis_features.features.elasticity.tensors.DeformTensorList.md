# matlantis_features.features.elasticity.tensors.DeformTensorList#

_class _matlantis_features.features.elasticity.tensors.DeformTensorList(_data : List[List[[DeformTensor](matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor "matlantis_features.features.elasticity.tensors.DeformTensor")]]_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#DeformTensorList)#
    

Bases: `object`

Tensor list object for the deformation tensors.

Methods

`__init__`(data) | Initialize an instance.  
---|---  
`from_strain_tensors`(strain_list) | Create deform tensor list from the strain tensor list.  
  
__init__(_data : List[List[[DeformTensor](matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor "matlantis_features.features.elasticity.tensors.DeformTensor")]]_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#DeformTensorList.__init__)#
    

Initialize an instance.

Parameters
    

**data** (_list_ _[__list_ _[_[_DeformTensor_](matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor "matlantis_features.features.elasticity.tensors.DeformTensor") _]__]_) – 2D-list of the deformation tensors.

_classmethod _from_strain_tensors(_strain_list : [StrainTensorList](matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList "matlantis_features.features.elasticity.tensors.StrainTensorList")_) → DeformTensorList[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#DeformTensorList.from_strain_tensors)#
    

Create deform tensor list from the strain tensor list.

Parameters
    

**strain_list** ([_StrainTensorList_](matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList "matlantis_features.features.elasticity.tensors.StrainTensorList")) – Strain tensor list object converted to the deform tensor list.

Returns
    

Deform tensor list object.

Return type
    

DeformTensorList

[ __ previous matlantis_features.features.elasticity.tensors.DeformTensor ](matlantis_features.features.elasticity.tensors.DeformTensor.html "previous page") [ next matlantis_features.features.elasticity.tensors.StrainTensor __](matlantis_features.features.elasticity.tensors.StrainTensor.html "next page")
