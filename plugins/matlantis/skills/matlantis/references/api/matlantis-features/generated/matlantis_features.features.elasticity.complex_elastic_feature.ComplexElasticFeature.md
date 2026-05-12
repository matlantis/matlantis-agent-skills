# matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature#

_class _matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature(_run_ase_opt_feature : [ASEOptFeature](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")_, _run_elastic_tensor_feature : [ElasticTensorFeature](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature "matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature")_, _post_elastic_feature : [PostElasticPropertiesFeature](matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature.html#matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature "matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature")_)[[source]](../_modules/matlantis_features/features/elasticity/complex_elastic_feature.html#ComplexElasticFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating elasticity related properties.

This feature calls both ElasticTensorFeature and PostElasticPropertiesFeature.

Methods

`__init__`(run_ase_opt_feature, ...) | Initialize an instance.  
---|---  
`__call__`(atoms) | Perform the feature's calculation.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_run_ase_opt_feature : [ASEOptFeature](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")_, _run_elastic_tensor_feature : [ElasticTensorFeature](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature "matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature")_, _post_elastic_feature : [PostElasticPropertiesFeature](matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature.html#matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature "matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature")_) → None[[source]](../_modules/matlantis_features/features/elasticity/complex_elastic_feature.html#ComplexElasticFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **run_ase_opt_feature** ([_ASEOptFeature_](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")) – Optimizer feature object used for elastic tensor calculation.

  * **run_elastic_tensor_feature** ([_ElasticTensorFeature_](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature "matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature")) – Feature object for calculating the elastic tensor.

  * **post_elastic_feature** ([_PostElasticPropertiesFeature_](matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature.html#matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature "matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature")) – Feature object for calculating the elastic properties.




__call__(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_) → [ComplexElasticFeatureResult](matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult.html#matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult "matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult")#
    

Perform the feature’s calculation.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – Atoms object containing the target system.

Returns
    

Dataclass object containing the calculation result.

Return type
    

[ComplexElasticFeatureResult](matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult.html#matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult "matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult")

attach_ctx(_ctx : Optional[Context] = None_) → None#
    

Attach the feature to matlantis_features.utils.Context.

Parameters
    

**ctx** (_Context_ _or_ _None_ _,__optional_) – The matlantis_features.utils.Context object. Defaults to None.

check_estimator_fn(_estimator_fn : Optional[Callable[[], Estimator]]_) → None#
    

Checks if the given estimator function is None and output a warning if so.

Parameters
    

**estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail. Defaults to None.

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

[ __ previous matlantis_features.features.elasticity.complex_elastic_feature ](matlantis_features.features.elasticity.complex_elastic_feature.html "previous page") [ next matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult __](matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult.html "next page")
