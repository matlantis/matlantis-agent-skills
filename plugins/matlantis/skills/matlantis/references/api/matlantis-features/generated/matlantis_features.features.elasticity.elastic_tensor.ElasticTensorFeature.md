# matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature#

_class _matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature(_diagonal_strains : List[float]_, _off_diagonal_strains : List[float]_, _optimizer : [ASEOptFeature](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")_, _check_linear_elems : Optional[List[Tuple[int, int]]] = None_, _zero_tol : float = 0.001_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _check_symmetry : bool = False_, _rtol : float = 0.01_, _atol : float = 1.0_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensorFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the elastic tensor for the specified atoms, by applying the specified strain.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/elastic_tensor.html).

Methods

`__init__`(diagonal_strains, ...[, ...]) | Initialize an instance.  
---|---  
`__call__`(matlantis_atoms) | Calculate the elastic tensor.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_elastic_tensor_symmetry`(atoms, ...) | Check correctness of elastic tensor according to the lattice symmetry.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`check_nonlinear`(elastic_tensor, ...) | Check non-linearity of stress-strain relation.  
`check_zero_elems`(elastic_tensor, ...) | Check if the elastic tensor elements are zero.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_diagonal_strains : List[float]_, _off_diagonal_strains : List[float]_, _optimizer : [ASEOptFeature](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")_, _check_linear_elems : Optional[List[Tuple[int, int]]] = None_, _zero_tol : float = 0.001_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _check_symmetry : bool = False_, _rtol : float = 0.01_, _atol : float = 1.0_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensorFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **diagonal_strains** (_list_ _[__float_ _]_) – Diagonal elements of the strain tensor.

  * **off_diagonal_strains** (_list_ _[__float_ _]_) – Off-diagonal elements of the strain tensor.

  * **optimizer** ([_ASEOptFeature_](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")) – Optimization method to optimize the system after applying the strain.

  * **check_linear_elems** (_list_ _[__tuple_ _[__int_ _,__int_ _]__] or_ _None_ _,__optional_) – Index list of the elastic tensor to check the linearlity. Default index list ((0, 0),(1, 1),(2, 2),(0, 1),(0, 2),(1, 2),(3, 3),(4, 4),(5, 5)) is used, if None is specified. Defaults to None.

  * **zero_tol** (_float_ _,__optional_) – Skip the linearity check if the elastic constant is smaller than this threshold. Defaults to 0.001.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **check_symmetry** (_bool_ _,__optional_) – Check the calculated elastic tensor is in accord with the type of lattice symmetry. For example, C11 = C22 and C14 = 0 in the cubic lattice. Please note that the two elements a and b are considered as equal if their difference is within a tolerance that us defined by rtol and atol, i.e. absolute(a - b) <= (atol + rtol * absolute(b)). If check_symmetry is enabled, the resultant elastic tensor corresponds to the rotated atomic coordinates after applying the symmetry operations, which might differ from the input atomic coordinates. Please call get_symmetry_lattice(atoms) if the rotated atomic coordinates are needed. In other words, if you do not want the atomic coordinates to be rotated from the original, please set check_symmetry to False. Default to False.

  * **rtol** (_float_ _,__optional_) – The relative tolerance of difference. The elastic tensor elements are regarded as the same if the relative difference is smaller than this value. Default to 0.01.

  * **atol** (_float_ _,__optional_) – The absolute tolerance of difference. The elastic tensor elements are regarded as the same if the absolute difference is smaller than this value. Default to 1.0

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.




__call__(_matlantis_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_) → [ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")#
    

Calculate the elastic tensor.

Parameters
    

**matlantis_atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – Atoms representing the system to calculate the elastic tensor.

Returns
    

Resulting elastic tensor object.

Return type
    

[ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")

attach_ctx(_ctx : Optional[Context] = None_) → None#
    

Attach the feature to matlantis_features.utils.Context.

Parameters
    

**ctx** (_Context_ _or_ _None_ _,__optional_) – The matlantis_features.utils.Context object. Defaults to None.

check_elastic_tensor_symmetry(_atoms : Atoms_, _elastic_tensor : [ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")_) → None[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensorFeature.check_elastic_tensor_symmetry)#
    

Check correctness of elastic tensor according to the lattice symmetry.

  1. Certain elements of elastic tensor should be zeros, except the triclinic lattice.

  2. Some elements of elastic tensor should be identical, e.g. C11 = C22 = C33 in cubic lattice.

  3. Some elements of elastic tensor should be opposite number, e.g. C16 = -C26 in tetragonal lattice.

  4. Certain relation exists between the elements, e.g. C66 = 1/2 (C11-C12) in hexagonal lattice.




Parameters
    

  * **atoms** (_ASEAtoms_) – Atoms representing the system to calculate the elastic tensor.

  * **elastic_tensor** ([_ElasticTensor_](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")) – Elastic tensor object containing the stress and strain data.




check_estimator_fn(_estimator_fn : Optional[Callable[[], Estimator]]_) → None#
    

Checks if the given estimator function is None and output a warning if so.

Parameters
    

**estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail. Defaults to None.

check_nonlinear(_elastic_tensor : [ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")_, _check_linear_elems : List[Tuple[int, int]]_) → None[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensorFeature.check_nonlinear)#
    

Check non-linearity of stress-strain relation.

If non linearity is detected, warning will be reported to the logger and self.report() will be called.

Parameters
    

  * **elastic_tensor** ([_ElasticTensor_](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")) – Elastic tensor object containing the stress and strain data.

  * **check_linear_elems** (_list_ _[__tuple_ _[__int_ _,__int_ _]__]_) – Index list of the tensor elements to check the linearity.




check_zero_elems(_elastic_tensor : [ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")_, _check_linear_elems : List[Tuple[int, int]]_, _zero_tol : float_) → None[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#ElasticTensorFeature.check_zero_elems)#
    

Check if the elastic tensor elements are zero.

If non zero emenets are detected, warning will be reported to the logger and self.report() will be called.

Parameters
    

  * **elastic_tensor** ([_ElasticTensor_](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")) – Elastic tensor object to be checked.

  * **check_linear_elems** (_list_ _[__tuple_ _[__int_ _,__int_ _]__]_) – Index list of the tensor elements to check the linearity. The elements specified in this list will not be checked, since these elements are expected to have non-zero elastic constants.

  * **zero_tol** (_float_) – Tolerance threshold for the zero check.




cost_estimate(_atoms : Optional[Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]] = None_) → FeatureCost#
    

Estimate the cost of the feature.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _None_ _,__optional_) – The input atoms. Defaults to None.

Returns
    

The cost of the feature.

Return type
    

FeatureCost

_classmethod _from_dict(_res : Dict[str, Any]_) → [FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")#
    

Construct a FeatureBase object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBase from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

[FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

get_savedir_from_ctx() → Path#
    

Get the temporary save directory from the context.

Returns
    

The temporary save directory .

Return type
    

pathlib.Path

init_scope() → Iterator[None]#
    

Context manager that enable to set attribution of the feature.

Returns
    

Init_scope context manager.

Return type
    

Iterator[None]

repeat(_n_repeat : int_) → Self#
    

Set the maximum number of times that allowed to run the __call__ function.

Parameters
    

**n_repeat** (_int_) – The maximum number of repeats.

Returns
    

The feature.

Return type
    

Self

to_dict() → Dict[str, Any]#
    

Dictionary representation of the FeatureBase.

Returns
    

A dict containing a serialized FeatureBase.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.elasticity.elastic_tensor ](matlantis_features.features.elasticity.elastic_tensor.html "previous page") [ next matlantis_features.features.elasticity.elastic_tensor.ElasticTensor __](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html "next page")
