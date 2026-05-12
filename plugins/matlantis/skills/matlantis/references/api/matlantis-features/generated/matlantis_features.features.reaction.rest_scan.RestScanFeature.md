# matlantis_features.features.reaction.rest_scan.RestScanFeature#

_class _matlantis_features.features.reaction.rest_scan.RestScanFeature(_scans : List[List[[ScanRestraint](matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint "restscan.restraints.ScanRestraint")]]_, _trajectory : Optional[str] = None_, _logfile : Optional[str] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 10_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

RestScan is a method for continuously changing molecular structures into different structures. This can be used to make initial estimates of reaction paths or to obtain alternative structures from existing molecular structures. It’s similar to the Gaussian scan function.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/rest_scan.html).

Methods

`__init__`(scans[, trajectory, logfile, ...]) | Initialize an instance.  
---|---  
`__call__`(atoms[, fmax, n_run, energy_stop, ...]) | Run the RestScan calculation.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_scans : List[List[[ScanRestraint](matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint "restscan.restraints.ScanRestraint")]]_, _trajectory : Optional[str] = None_, _logfile : Optional[str] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _logger_interval : int = 10_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **scans** (_list_ _[__list_ _[_[_ScanRestraint_](matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint "matlantis_features.features.reaction.rest_scan.ScanRestraint") _]__]_) – The scans is given as a list of lists of scans. The number of scans to be applied is increased in the order of the outer list.

  * **trajectory** (_str_ _or_ _None_ _,__optional_) – Attach trajectory object. If _trajectory_ is a string a Trajectory will be constructed. Use _None_ for no trajectory. Defaults to None.

  * **logfile** (_str_ _or_ _None_ _,__optional_) – If logfile is a string, a file with that name will be opened. Use ‘-’ for stdout. Defaults to None.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **logger_interval** (_int_ _,__optional_) – The interval of when to print out logger information. Defaults to 10.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.




__call__(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _fmax : float = 0.01_, _n_run : int = 10000_, _energy_stop : Optional[float] = None_, _local_minima_de : float = 0.01_, _local_minima_n : int = 1_, _attach_methods : Optional[List[Callable[[], None]]] = None_) → [RestScanFeatureResult](matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult "matlantis_features.features.reaction.rest_scan.RestScanFeatureResult")#
    

Run the RestScan calculation.

Parameters
    

  * **atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – The input structure.

  * **fmax** (_float_ _,__optional_) – The maximum force (in eV/Angstrom) for convergence of scan. Defaults to 0.01.

  * **n_run** (_int_ _,__optional_) – The maximum number of scan steps performed to try to reach force convergence (specified by fmax). Defaults to 10000.

  * **energy_stop** (_float_ _or_ _None_ _,__optional_) – Exit when the potential energy of the atoms is higher than this value. Defaults to None.

  * **local_minima_de** (_float_ _,__optional_) – Local minimas that are lower than local maxima by this value are counted as local minimas. Defaults to 0.01.

  * **local_minima_n** (_int_ _,__optional_) – Exit when the trajectory reach n’th local minima. Defaults to 1.

  * **attach_methods** (_list_ _[__- > None_ _] or_ _None_ _,__optional_) – Functions to be called by the restscan. A list of callable functions. If None, no functions will be attached. Defaults to None.



Returns
    

The calculation result.

Return type
    

[RestScanFeatureResult](matlantis_features.features.reaction.rest_scan.RestScanFeatureResult.html#matlantis_features.features.reaction.rest_scan.RestScanFeatureResult "matlantis_features.features.reaction.rest_scan.RestScanFeatureResult")

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

[ __ previous RestScan Features ](../matlantis_features.features.reaction.rest_scan.html "previous page") [ next matlantis_features.features.reaction.rest_scan.ScanDihedral __](matlantis_features.features.reaction.rest_scan.ScanDihedral.html "next page")
