# matlantis_features.features.elasticity.tensors.StrainTensor#

_class _matlantis_features.features.elasticity.tensors.StrainTensor(_vix : int_, _strain : float_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensor)#
    

Bases: `_Tensor3x3`

Strain tensor object.

Methods

`__init__`(vix, strain) | Initialize an instance.  
---|---  
`get_by_voigt`(vix) | Get tensor element by voigt index.  
`set_by_voigt`(vix, value) | Set tensor element by voigt index.  
`to_3x3_ndarray`() | Get 3x3 matrix representation of this tensor.  
`to_voigt_ndarray`() | Get 6x1 Voigt representation of this tensor.  
  
__init__(_vix : int_, _strain : float_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StrainTensor.__init__)#
    

Initialize an instance.

Parameters
    

  * **vix** (_int_) – Voigt index.

  * **strain** (_float_) – Strain value specified by the voigt index, vix.




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

[ __ previous matlantis_features.features.elasticity.tensors.DeformTensorList ](matlantis_features.features.elasticity.tensors.DeformTensorList.html "previous page") [ next matlantis_features.features.elasticity.tensors.StrainTensorList __](matlantis_features.features.elasticity.tensors.StrainTensorList.html "next page")
