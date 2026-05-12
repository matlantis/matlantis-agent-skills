# matlantis_features.features.common.vibration.VibrationFeature#

_class _matlantis_features.features.common.vibration.VibrationFeature(_delta : float = 0.01_, _nfree : int = 2_, _method : str = 'standard'_, _direction : str = 'central'_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _num_threads : int = 20_)[[source]](../_modules/matlantis_features/features/common/vibration.html#VibrationFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating molecular vibration.

The vibration properties, including the force constant, vibration frequencies and vibration modes are calculated with the finite displacement method.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/vibration.html).

Methods

`__init__`([delta, nfree, method, direction, ...]) | Initialize an instance.  
---|---  
`__call__`(atoms[, indices]) | Run the vibration calculation.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_delta : float = 0.01_, _nfree : int = 2_, _method : str = 'standard'_, _direction : str = 'central'_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _num_threads : int = 20_) → None[[source]](../_modules/matlantis_features/features/common/vibration.html#VibrationFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **delta** (_float_ _,__optional_) – The distance of finite displacement applied to each atom. The unit is Angstrom. Defaults to 0.01.

  * **nfree** (_int_ _,__optional_) – The number of displacements per atom. Must be 2 or 4. Defaults to 2.

  * **method** (_str_ _,__optional_) – The method to calculate the force constant from the displaced supercells and their forces. This parameter can be ‘Standard’ or ‘Frederiksen’. Defaults to ‘standard’.

  * **direction** (_str_ _,__optional_) – This parameter discribe how the displacment is found. Must be When ‘direction=forward’/’direction=backward’, each atom is displaced along the positive/negative direction of the axis. When ‘direction=center’, each atom is displaced along both positive and negative directions of the axis. Defaults to ‘central’.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.

  * **num_threads** (_int_ _,__optional_) – Number of threads for parallel calculation. Defaults to 10.




__call__(_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_, _indices : Optional[List[int]] = None_) → [VibrationFeatureResult](matlantis_features.features.common.vibration.VibrationFeatureResult.html#matlantis_features.features.common.vibration.VibrationFeatureResult "matlantis_features.features.common.vibration.VibrationFeatureResult")#
    

Run the vibration calculation.

Parameters
    

  * **atoms** ([_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_) – The input structure.

  * **indices** (_list_ _[__int_ _] or_ _None_ _,__optional_) – List of indices of atoms allowed to vibrate. If ‘None’ is provided, all atoms are allowed to vibrate. Defaults to None.



Returns
    

The vibration calculation results, which include the force constant, vibration frequencies and vibration modes.

Return type
    

[VibrationFeatureResult](matlantis_features.features.common.vibration.VibrationFeatureResult.html#matlantis_features.features.common.vibration.VibrationFeatureResult "matlantis_features.features.common.vibration.VibrationFeatureResult")

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

[ __ previous matlantis_features.features.common.vibration ](matlantis_features.features.common.vibration.html "previous page") [ next matlantis_features.features.common.vibration.VibrationFeatureResult __](matlantis_features.features.common.vibration.VibrationFeatureResult.html "next page")
