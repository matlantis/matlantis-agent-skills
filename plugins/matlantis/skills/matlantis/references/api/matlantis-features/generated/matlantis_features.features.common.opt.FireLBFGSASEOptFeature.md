# matlantis_features.features.common.opt.FireLBFGSASEOptFeature#

_class _matlantis_features.features.common.opt.FireLBFGSASEOptFeature(_n_run : int = 200_, _fmax : float = 0.01_, _switch : float = 0.05_, _switch_decrease_rate : float = 0.9_, _switch_patience_fire : int = 10_, _switch_patience_lbfgs : int = 10_, _attach_methods : Optional[List[Tuple[Callable[[], None], int]]] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 10_, _filter : Union[bool, [MatlantisFilter](../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter "matlantis_features.filters.MatlantisFilter")] = False_, _trajectory : Optional[str] = None_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _maxstep_fire : float = 0.2_, _maxstep_lbfgs : float = 0.2_)[[source]](../_modules/matlantis_features/features/common/opt.html#FireLBFGSASEOptFeature)#
    

Bases: [`ASEOptFeature`](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")

The matlantis-feature that combines FIRE and LBFGS optimization algorithms.

Methods

`__init__`([n_run, fmax, switch, ...]) | Initialize an instance.  
---|---  
`__call__`(atoms) | Call function.  
`attach`(extn, interval) | Used for attaching a function to the ASE optimizer.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FireLBFGSASEOptFeature.  
  
__init__(_n_run : int = 200_, _fmax : float = 0.01_, _switch : float = 0.05_, _switch_decrease_rate : float = 0.9_, _switch_patience_fire : int = 10_, _switch_patience_lbfgs : int = 10_, _attach_methods : Optional[List[Tuple[Callable[[], None], int]]] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 10_, _filter : Union[bool, [MatlantisFilter](../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter "matlantis_features.filters.MatlantisFilter")] = False_, _trajectory : Optional[str] = None_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _maxstep_fire : float = 0.2_, _maxstep_lbfgs : float = 0.2_)[[source]](../_modules/matlantis_features/features/common/opt.html#FireLBFGSASEOptFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **n_run** (_int_ _,__optional_) – The maximum number of optimization steps performed to try to reach force convergence (specified by fmax). Defaults to 200.

  * **fmax** (_float_ _,__optional_) – The maximum force (in eV/Angstrom) for convergence of optimization. Defaults to 0.01.

  * **attach_methods** (_list_ _[__tuple_ _[__- > None_ _,__int_ _]__] or_ _None_ _,__optional_) – Functions to be called by the optimizer. A list of tuples must be specified, where the first element of each tuple is a callable function, and the second element is an integer N. The optimizer will call the function every N steps. If None, no functions will be attached. Defaults to None.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **logger_interval** (_int_ _,__optional_) – The interval of when to print out logger information. Defaults to 10.

  * **filter** (_bool_ _or_[ _MatlantisFilter_](../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter "matlantis_features.filters.MatlantisFilter") _,__optional_) – The MatlantisFilter used in the optimization. If True, the default filter, i.e. UnitCellASEFilter, will be used for the cell shape optimization. If False, no filter is used. Note: The filter is only applicable in the optimization of ASEAtoms and MatlantisAtoms, not NEB class. Defaults to False.

  * **trajectory** (_str_ _or_ _None_ _,__optional_) – Filename for trajectory file. If None, no trajectory file will be saved. Defaults to None.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.

  * **maxstep** (_float_ _,__optional_) – The maximum distance an atom can move per iteration. The unit is angstrom. Defaults to 0.2.




__call__(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), NEB, Filter]_) → [OptFeatureResult](matlantis_features.features.common.opt.OptFeatureResult.html#matlantis_features.features.common.opt.OptFeatureResult "matlantis_features.features.common.opt.OptFeatureResult")#
    

Call function.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _NEB_ _or_ _Filter_) – The system to be optimized. In addition to ASEAtoms and MatlantisAtoms, NEB (Matlantis) and Filter (ASE) objects are also supported. Note: the support of ASE Filter class will be deprecated in the future. Please use the ‘filter’ parameter in the __init__ function.

Returns
    

A dataclass containing information about the result of optimization.

Return type
    

[OptFeatureResult](matlantis_features.features.common.opt.OptFeatureResult.html#matlantis_features.features.common.opt.OptFeatureResult "matlantis_features.features.common.opt.OptFeatureResult")

attach(_extn : Callable[[], None]_, _interval : int_) → None#
    

Used for attaching a function to the ASE optimizer.

Parameters
    

  * **extn** (_- > None_) – The function to be attached.

  * **interval** (_int_) – The interval (in optimization steps) of when the optimizer will call the function.




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

to_dict() → Dict[str, Any][[source]](../_modules/matlantis_features/features/common/opt.html#FireLBFGSASEOptFeature.to_dict)#
    

Dictionary representation of the FireLBFGSASEOptFeature.

Returns
    

A dict containing a serialized FireLBFGSASEOptFeature.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.common.opt.FireASEOptFeature ](matlantis_features.features.common.opt.FireASEOptFeature.html "previous page") [ next matlantis_features.features.common.opt.LBFGSASEOptFeature __](matlantis_features.features.common.opt.LBFGSASEOptFeature.html "next page")
