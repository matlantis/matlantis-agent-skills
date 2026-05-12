# matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature#

_class _matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature[[source]](../_modules/matlantis_features/features/md/post_features/emd_viscosity.html#PostEMDViscosityFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the viscosity and related properties, using the equilibrium MD trajectory.

References

This implementationis based on the literature: J. Chem. Phys. 106, 9327 (1997); <https://aip.scitation.org/doi/10.1063/1.474002>

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/viscosity_emd.html).

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(md_results[, init_time, ...]) | Calculate the viscosity from the MD trajectory.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None#
    

Initialize an instance.

__call__(_md_results : [MDFeatureResult](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")_, _init_time : float = 0.0_, _number_of_segments : int = 1_) → [PostEMDViscosityFeatureResult](matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult "matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult")#
    

Calculate the viscosity from the MD trajectory.

Parameters
    

  * **md_results** ([_MDFeatureResult_](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")) – MD feature’s result object containing the trajectory.

  * **init_time** (_float_ _,__optional_) – Trajectory frames before this time is ignored for the calculation. Defaults to 0.0 (i.e., all frames are used for calculation).

  * **number_of_segments** (_int_ _,__optional_) – Divides the given trajectory in to segments to allow statistical analysis. Defaults to 1.



Returns
    

Result dataclass object.

Return type
    

[PostEMDViscosityFeatureResult](matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult "matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult")

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

[ __ previous matlantis_features.features.md.post_features.emd_viscosity ](matlantis_features.features.md.post_features.emd_viscosity.html "previous page") [ next matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult __](matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html "next page")
