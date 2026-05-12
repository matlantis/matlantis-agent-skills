# matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature#

_class _matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature[[source]](../_modules/matlantis_features/features/md/post_features/thermal_expansion.html#PostMDThermalExpansionFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the thermal expansion from the MD trajectories at different temperatures.

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(md_results_list, temperatures[, ...]) | Call function.  
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

__call__(_md_results_list : List[[MDFeatureResult](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")]_, _temperatures : List[float]_, _init_time : float = 0.0_, _fitting_order : int = 3_) → [PostMDThermalExpansionFeatureResult](matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.html#matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult "matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult")#
    

Call function.

Parameters
    

  * **md_results_list** (_list_ _[_[_MDFeatureResult_](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult") _]_) – MD feature’s result object containing the trajectory.

  * **temperatures** (_list_ _[__float_ _]_) – List of the temperatures the simulations were performed.

  * **init_time** (_float_ _,__optional_) – Trajectory frames before this time is ignored for the calculation. Defaults to 0.0 (i.e., all frames are used for calculation).

  * **fitting_order** (_int_ _,__optional_) – The order used for the polynomial fitting of the T-H relationship. Defaults to 3.



Returns
    

Result dataclass object.

Return type
    

[PostMDThermalExpansionFeatureResult](matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.html#matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult "matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult")

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

[ __ previous matlantis_features.features.md.post_features.thermal_expansion ](matlantis_features.features.md.post_features.thermal_expansion.html "previous page") [ next matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult __](matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.html "next page")
