# matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature#

_class _matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#PostPhononQHAGibbsFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating Gibbs free energy from Quasi-harmonic approximation.

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(force_constant_list, kpts, ...) | Calculates the Gibbs free energy from a list of force constants.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#PostPhononQHAGibbsFeature.__init__)#
    

Initialize an instance.

__call__(_force_constant_list : List[[ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")]_, _kpts : List[int]_, _temperatures : List[float]_, _pressures : List[float]_) → [PostPhononQHAGibbsFeatureResult](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult "matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult")#
    

Calculates the Gibbs free energy from a list of force constants.

Parameters
    

  * **force_constant_list** (_list_ _[_[_ForceConstantFeatureResult_](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult") _]_) – A list of force constant. The force constants are obtained by the ForceConstantFeature. The unit cell should be same except the volume is different. The len(force_constant_list) must be larger than 3.

  * **kpts** (_list_ _[__int_ _]_) – The number of mesh points along each axis in the Helmholtz free energy calculation

  * **temperatures** (_list_ _[__float_ _]_) – The temperatures at which the Gibbs free energy will be calculated.

  * **pressures** (_list_ _[__float_ _]_) – The pressures at which the Gibbs free energy will be calculated.



Returns
    

The calculation results of Gibbs free energy calculation.

Return type
    

[PostPhononQHAGibbsFeatureResult](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult "matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult")

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

[ __ previous matlantis_features.features.phonon.qha_gibbs_free_energy ](matlantis_features.features.phonon.qha_gibbs_free_energy.html "previous page") [ next matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult __](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html "next page")
