# matlantis_features.features.elasticity.elastic_tensor.ElasticTensor#

_class _matlantis_features.features.elasticity.elastic_tensor.ElasticTensor(_stress_data : ndarray_, _strain_data : ndarray_)[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensor)#
    

Bases: `object`

Elastic tensor object calculated from the stress-strain relation.

Methods

`__init__`(stress_data, strain_data) | Initialize an instance.  
---|---  
`fitting_result_linear`(i, j[, threshold]) | Check the linear relation between stress and strain values.  
`from_stress_strain`(stress, strain) | Create elastic tensor object from stress and strain tensor objects.  
  
__init__(_stress_data : ndarray_, _strain_data : ndarray_) → None[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensor.__init__)#
    

Initialize an instance.

Parameters
    

  * **stress_data** (_np.ndarray_) – Stress tensor for the calculation.

  * **strain_data** (_np.ndarray_) – Strain tensor for the calculation.




fitting_result_linear(_i : int_, _j : int_, _threshold : float = 0.5_) → bool[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensor.fitting_result_linear)#
    

Check the linear relation between stress and strain values.

Parameters
    

  * **i** (_int_) – i-index for the tensor.

  * **j** (_int_) – j-index for the tensor.

  * **threshold** (_float_ _,__optional_) – Linearlity criterion threshold. Defaults to 0.5.



Returns
    

Returns True if the relation satisfy the linearlity criterion threshold.

Return type
    

bool

_classmethod _from_stress_strain(_stress : [StressTensorList](matlantis_features.features.elasticity.tensors.StressTensorList.html#matlantis_features.features.elasticity.tensors.StressTensorList "matlantis_features.features.elasticity.tensors.StressTensorList")_, _strain : [StrainTensorList](matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList "matlantis_features.features.elasticity.tensors.StrainTensorList")_) → ElasticTensor[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensor.from_stress_strain)#
    

Create elastic tensor object from stress and strain tensor objects.

Parameters
    

  * **stress** ([_StressTensorList_](matlantis_features.features.elasticity.tensors.StressTensorList.html#matlantis_features.features.elasticity.tensors.StressTensorList "matlantis_features.features.elasticity.tensors.StressTensorList")) – Stress tensor list object.

  * **strain** ([_StrainTensorList_](matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList "matlantis_features.features.elasticity.tensors.StrainTensorList")) – Strain tensor list object.



Returns
    

Resulting elastic tensor object.

Return type
    

ElasticTensor

[ __ previous matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature ](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html "previous page") [ next matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms __](matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms.html "next page")
