# matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature#

_class _matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature(_integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _n_run : int_, _rnemd_n_slab : int = 20_, _rnemd_interval : int = 100_, _init_time : float = 0.0_, _md_feature_args : Optional[Dict[str, Any]] = None_)[[source]](../_modules/matlantis_features/features/md/post_features/nemd_viscosity.html#ComplexNEMDViscosityFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating viscosity related properties.

This feature calls both MDFeature with RNEMDExtension and PostNEMDViscosityFeature.

Methods

`__init__`(integrator, n_run[, rnemd_n_slab, ...]) | Initialize an instance.  
---|---  
`__call__`(system) | Run the rNEMD simulation and calculate the viscosity.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _n_run : int_, _rnemd_n_slab : int = 20_, _rnemd_interval : int = 100_, _init_time : float = 0.0_, _md_feature_args : Optional[Dict[str, Any]] = None_)[[source]](../_modules/matlantis_features/features/md/post_features/nemd_viscosity.html#ComplexNEMDViscosityFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation run.

  * **n_run** (_int_) – Number of time steps of the MD simulation.

  * **rnemd_n_slab** (_int_ _,__optional_) – Number of the slabs used for the rNEMD calculation. Defaults to 20.

  * **rnemd_interval** (_int_ _,__optional_) – Number of the timestep interval for the rNEMD calculation. Defaults to 100.

  * **init_time** (_float_ _,__optional_) – Trajectory frames before this time (in fs) is ignored for the calculation. Defaults to 0.0 (i.e., all frames are used for calculation).

  * **md_feature_args** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Additional options for MDFeature




__call__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_) → [ComplexNEMDViscosityFeatureResult](matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult "matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult")#
    

Run the rNEMD simulation and calculate the viscosity.

Parameters
    

**system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system to calculate the viscosity.

Returns
    

Result dataclass object.

Return type
    

[ComplexNEMDViscosityFeatureResult](matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult "matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult")

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

[ __ previous matlantis_features.features.md.post_features.nemd_viscosity ](matlantis_features.features.md.post_features.nemd_viscosity.html "previous page") [ next matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature __](matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature.html "next page")
