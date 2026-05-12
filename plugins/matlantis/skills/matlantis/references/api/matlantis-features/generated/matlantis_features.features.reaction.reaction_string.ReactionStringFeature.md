# matlantis_features.features.reaction.reaction_string.ReactionStringFeature#

_class _matlantis_features.features.reaction.reaction_string.ReactionStringFeature(_optimize_is : bool = True_, _optimize_fs : bool = True_, _fmax_rp : float = 0.05_, _fmax_ts : float = 0.05_, _fmax_rd : float = 0.01_, _fmax_eq : float = 0.001_, _local_extrema_tol : float = 0.001_, _interesting : Optional[Callable[[List[List[Atoms]]], List[bool]]] = None_, _timeout : Optional[timedelta] = None_, _logfile_rp : Optional[str] = None_, _logfile_eq : Optional[str] = None_, _observers_rp : Optional[List[Callable[[], None]]] = None_, _observers_eq : Optional[List[Callable[[], None]]] = None_, _remove_rotation_and_translation : Optional[bool] = None_, _max_workers : Optional[int] = None_, _optimizer_eq : Union[OptimizerBuilder, str] = 'FIRELBFGS'_, _optimizer_rp : Union[OptimizerBuilder, str] = 'FIRES'_, _k : float = 0.1_, _dx_ts : float = 0.02_, _dx_rp : float = 0.3_, _dx_rm : float = 1.0_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _dump_directory : Optional[Path] = None_, _dump_frequency : Optional[timedelta] = None_)[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

Reliable string method that surely connect IS and FS.

All local minima are true local minima. All local maxima is the result of dimer method.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/reaction_string.html).

Methods

`__init__`([optimize_is, optimize_fs, ...]) | Initialize an instance.  
---|---  
`__call__`(images) | Get the reaction path between the initial and end images.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the ReactionStringFeature.  
  
__init__(_optimize_is : bool = True_, _optimize_fs : bool = True_, _fmax_rp : float = 0.05_, _fmax_ts : float = 0.05_, _fmax_rd : float = 0.01_, _fmax_eq : float = 0.001_, _local_extrema_tol : float = 0.001_, _interesting : Optional[Callable[[List[List[Atoms]]], List[bool]]] = None_, _timeout : Optional[timedelta] = None_, _logfile_rp : Optional[str] = None_, _logfile_eq : Optional[str] = None_, _observers_rp : Optional[List[Callable[[], None]]] = None_, _observers_eq : Optional[List[Callable[[], None]]] = None_, _remove_rotation_and_translation : Optional[bool] = None_, _max_workers : Optional[int] = None_, _optimizer_eq : Union[OptimizerBuilder, str] = 'FIRELBFGS'_, _optimizer_rp : Union[OptimizerBuilder, str] = 'FIRES'_, _k : float = 0.1_, _dx_ts : float = 0.02_, _dx_rp : float = 0.3_, _dx_rm : float = 1.0_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _dump_directory : Optional[Path] = None_, _dump_frequency : Optional[timedelta] = None_) → None[[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **optimize_is** (_bool_ _,__optional_) – Optimize initial structure. Defaults to True.

  * **optimize_fs** (_bool_ _,__optional_) – Optimize final structure. Defaults to True.

  * **fmax_rp** (_float_ _,__optional_) – Convergence criteria for reaction path. Defaults to 0.05.

  * **fmax_ts** (_float_ _,__optional_) – Convergence criteria for transition states. If math.isinf(fmax_ts), climbing algorithm will not be executed while optimizing. Defaults to 0.05.

  * **fmax_rd** (_float_ _,__optional_) – Convergence criterial for rate determining step (highest TS). If math.isinf(fmax_rd), climbing algorithm will not be executed while optimizing. Defaults to 0.01.

  * **fmax_eq** (_float_ _,__optional_) – Convergence criteria for local minimas. Defaults to 0.001.

  * **local_extrema_tol** (_float_ _,__optional_) – Local extrema smaller than this will be ignored. Defaults to 1e-3.

  * **interesting** (_list_ _[__list_ _[__Atoms_ _]__]__- > list_ _[__bool_ _] or_ _None_ _,__optional_) – Interesting path or not. Defaults to None.

  * **timeout** (_timedelta_ _or_ _None_ _,__optional_) – Timeout. Defaults to None.

  * **logfile_rp** (_str_ _or_ _None_ _,__optional_) – Logfile for reaction path finding algorithms. If _logfile_ is a string, a file with that name will be opened. Use ‘-’ for stdout. Defaults to None.

  * **logfile_eq** (_str_ _or_ _None_ _,__optional_) – Logfile for local minima optimization. If _logfile_ is a string, a file with that name will be opened. Use ‘-’ for stdout. Defaults to None.

  * **observers_rp** (_list_ _[__- > None_ _] or_ _None_ _,__optional_) – Called in each step of reaction path optimzation. Defaults to None.

  * **observers_eq** (_list_ _[__- > None_ _] or_ _None_ _,__optional_) – Called after local minimas are optimized. Defaults to None.

  * **remove_rotation_and_translation** (_bool_ _or_ _None_ _,__optional_) – True actives NEB-TR for removing translation and rotation during NEB. By default applied non-periodic systems Defaults to None.

  * **max_workers** (_int_ _or_ _None_ _,__optional_) – Max number of concurrent threads. If it is not specified, it is automatically determined based on the number of input atoms.

  * **optimizer_eq** (_OptimizerBuilder_ _or_ _str_ _,__optional_) – Optimizer generator for structural optimization. Defaults to “FIRELBFGS”. Passing OptimizerBuilder type objects will be deprecated in a future version.

  * **optimizer_rp** (_OptimizerBuilder_ _or_ _str_ _,__optional_) – Optimizer generator for reaction path optimization. Defaults to “FIRE”. Passing OptimizerBuilder type objects will be deprecated in a future version.

  * **k** (_float_ _,__optional_) – NEB spring constant for String and NEB hybrid algorithm. Defaults to 0.1.

  * **dx_ts** (_float_ _,__optional_) – Distance between images for climbing string. Defaults to 0.02.

  * **dx_rp** (_float_ _,__optional_) – Distance between images for reaction path. Defaults to 0.3.

  * **dx_rm** (_float_ _,__optional_) – Images closer than this value will be linearly interpolated to skip connection. Defaults to 1.0.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail. Defaults to None.

  * **dump_directory** (_Path_ _or_ _None_ _,__optional_) – Folder to save images during optimization. If dump_directory is not None, you must specify dump_frequency as well.

  * **dump_frequency** (_timedelta_ _or_ _None_ _,__optional_) – Interval to dump images on dump_directory. If dump_frequency is not None, you must specify dump_directory as well.




__call__(_images : List[List[Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]]]_) → [ReactionStringFeatureResult](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult "matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult")#
    

Get the reaction path between the initial and end images.

Parameters
    

**images** (_list_ _[__list_ _[_[_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_ _]__]_) – Images defining path from initial to final state.

Returns
    

The calculation result.

Return type
    

[ReactionStringFeatureResult](matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult.html#matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult "matlantis_features.features.reaction.reaction_string.ReactionStringFeatureResult")

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

to_dict() → Dict[str, Any][[source]](../_modules/matlantis_features/features/reaction/reaction_string.html#ReactionStringFeature.to_dict)#
    

Dictionary representation of the ReactionStringFeature.

Returns
    

A dict containing a serialized ReactionStringFeature.

Return type
    

dict[str, Any]

[ __ previous ReactionString Features ](../matlantis_features.features.reaction.reaction_string.html "previous page") [ next matlantis_features.features.reaction.reaction_string __](matlantis_features.features.reaction.reaction_string.html "next page")
