# matlantis_features.features.md.md.MDFeature#

_class _matlantis_features.features.md.md.MDFeature(_integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _n_run : int_, _checkpoint_file_name : Optional[str] = None_, _checkpoint_freq : Optional[int] = 10_, _traj_file_name : Optional[str] = None_, _traj_freq : Optional[int] = 1_, _traj_props : Optional[List[str]] = None_, _traj_append : bool = False_, _traj_init_frame : bool = False_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 100_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/md/md.html#MDFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for molecular dynamics simulation.

Methods

`__init__`(integrator, n_run[, ...]) | Initialize an instance.  
---|---  
`__call__`(system[, extensions]) | Run the MD simulation.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`get_temp_dir`() | Create temperory directory.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _n_run : int_, _checkpoint_file_name : Optional[str] = None_, _checkpoint_freq : Optional[int] = 10_, _traj_file_name : Optional[str] = None_, _traj_freq : Optional[int] = 1_, _traj_props : Optional[List[str]] = None_, _traj_append : bool = False_, _traj_init_frame : bool = False_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 100_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/md/md.html#MDFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator object used for the simulation.

  * **n_run** (_int_) – Number or the steps to be run.

  * **checkpoint_file_name** (_str_ _or_ _None_ _,__optional_) – Checkpoint file name. If this value is None and checkpoint_freq>0, checkpoint file name is automatically generated. If relative path name is given, the file is stored in the current working directory. Defaults to None.

  * **checkpoint_freq** (_int_ _or_ _None_ _,__optional_) – Checkpoint creation frequency. If this value is None or non positive, no checkpoit file is created. Defaults to 10.

  * **traj_file_name** (_str_ _or_ _None_ _,__optional_) – Trajectory file name. If this value is None and traj_freq>0, trajectory file name is automatically generated. If relative path name is given, the file is stored in the current working directory. Defaults to None.

  * **traj_freq** (_int_ _or_ _None_ _,__optional_) – Trajectory saving frequency. If this value is None or non positive, no trajectory file is created. Defaults to 1.

  * **traj_props** (_list_ _[__str_ _] or_ _None_ _,__optional_) – Names of the properties stored in the trajectory file. If None, energy, forces, and stress properties are stored. Defaults to None.

  * **traj_append** (_bool_ _,__optional_) – If True and trajectory file exists, trajectory will be appended to the existing file. Defaults to False.

  * **traj_init_frame** (_bool_ _,__optional_) – If True, the initial structure of the simulation is also stored in the trajectory file. Defaults to False.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **logger_interval** (_int_ _,__optional_) – The interval of when to print out logger information. Defaults to 100.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.




__call__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _extensions : Optional[List[Tuple[[MDExtensionBase](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase"), int]]] = None_) → [MDFeatureResult](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")#
    

Run the MD simulation.

Parameters
    

  * **system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system of the MD simulation.

  * **extensions** (_MDExtensionList_ _or_ _None_ _,__optional_) – Extension list (tuple of extension object and interval) used for the simulation. Defaults to None.



Returns
    

Result object for the MD simulation.

Return type
    

[MDFeatureResult](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")

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

get_temp_dir() → Path[[source]](../_modules/matlantis_features/features/md/md.html#MDFeature.get_temp_dir)#
    

Create temperory directory.

Returns
    

The location of temperory directory.

Return type
    

Path

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

[ __ previous matlantis_features.features.md.md ](matlantis_features.features.md.md.html "previous page") [ next matlantis_features.features.md.md.MDFeatureResult __](matlantis_features.features.md.md.MDFeatureResult.html "next page")
