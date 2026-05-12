# matlantis_features.features.elasticity.tensors.StressTensor#

_class _matlantis_features.features.elasticity.tensors.StressTensor(_stress : ndarray_)[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StressTensor)#
    

Bases: `_TensorVoigt`

Stress tensor object.

Methods

`__init__`(stress) | Initialize an instance.  
---|---  
`get_by_voigt`(vix) | Get tensor element by voigt index.  
`set_by_voigt`(vix, value) | Set tensor element by voigt index.  
`to_3x3_ndarray`() | Get 3x3 matrix representation of this tensor.  
`to_voigt_ndarray`() | Get 6x1 Voigt representation of this tensor.  
  
__init__(_stress : ndarray_) → None[[source]](../_modules/matlantis_features/features/elasticity/tensors.html#StressTensor.__init__)#
    

Initialize an instance.

Parameters
    

**stress** (_np.ndarray_) – Stress tensor value in Voigt’s notation.

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

[ __ previous matlantis_features.features.elasticity.tensors.StrainTensorList ](matlantis_features.features.elasticity.tensors.StrainTensorList.html "previous page") [ next matlantis_features.features.elasticity.tensors.StressTensorList __](matlantis_features.features.elasticity.tensors.StressTensorList.html "next page")
