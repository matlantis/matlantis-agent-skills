# matlantis_features.features.elasticity.tensors.DeformTensor#

_class _matlantis_features.features.elasticity.tensors.DeformTensor(_strain_matrix : [StrainTensor](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor")_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#DeformTensor)#
    

Bases: `_Tensor3x3`

Deformation tensor object.

Methods

`__init__`(strain_matrix) | Initialize an instance.  
---|---  
`get_by_voigt`(vix) | Get tensor element by voigt index.  
`set_by_voigt`(vix, value) | Set tensor element by voigt index.  
`to_3x3_ndarray`() | Get 3x3 matrix representation of this tensor.  
`to_voigt_ndarray`() | Get 6x1 Voigt representation of this tensor.  
  
__init__(_strain_matrix : [StrainTensor](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor")_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#DeformTensor.__init__)#
    

Initialize an instance.

Parameters
    

**strain_matrix** ([_StrainTensor_](matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor "matlantis_features.features.elasticity.tensors.StrainTensor")) – Strain tensor object converted to the deformation tensor.

get_by_voigt(_vix : int_) → float#
    

Get tensor element by voigt index.

Parameters
    

**vix** (_int_) – Voigt index (0-5).

Returns
    

Value of the tensor of the specified index.

Return type
    

float

set_by_voigt(_vix : int_, _value : float_) → None#
    

Set tensor element by voigt index.

Parameters
    

  * **vix** (_int_) – Voigt index (0-5).

  * **value** (_float_) – Value to be stored in the tensor.




to_3x3_ndarray() → ndarray#
    

Get 3x3 matrix representation of this tensor.

Returns
    

numpy array of (3,3) shape.

Return type
    

np.ndarray

to_voigt_ndarray() → ndarray#
    

Get 6x1 Voigt representation of this tensor.

Returns
    

numpy array of (6,1) shape.

Return type
    

np.ndarray

[ __ previous matlantis_features.features.elasticity.tensors ](matlantis_features.features.elasticity.tensors.html "previous page") [ next matlantis_features.features.elasticity.tensors.DeformTensorList __](matlantis_features.features.elasticity.tensors.DeformTensorList.html "next page")
