# matlantis_features.features.common.vibration.FastVibrations#

_class _matlantis_features.features.common.vibration.FastVibrations(_delta : float = 0.01_, _nfree : int = 2_, _method : str = 'standard'_, _direction : str = 'central'_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _num_threads : int = 10_)[[source]](../_modules/matlantis_features/features/common/vibration.html#FastVibrations)#
    

Bases: `object`

The class for calculating the force constant, vibration frequency and vibration mode using the finite displacement method.

Methods

`__init__`([delta, nfree, method, direction, ...]) | Initialize an instance.  
---|---  
`run`(atoms[, indices]) | Run vibration calculation.  
  
__init__(_delta : float = 0.01_, _nfree : int = 2_, _method : str = 'standard'_, _direction : str = 'central'_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_, _num_threads : int = 10_) → None[[source]](../_modules/matlantis_features/features/common/vibration.html#FastVibrations.__init__)#
    

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




run(_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_, _indices : Optional[List[int]] = None_) → [FastVibrationsResult](matlantis_features.features.common.vibration.FastVibrationsResult.html#matlantis_features.features.common.vibration.FastVibrationsResult "matlantis_features.features.common.vibration.FastVibrationsResult")[[source]](../_modules/matlantis_features/features/common/vibration.html#FastVibrations.run)#
    

Run vibration calculation.

Parameters
    

  * **atoms** ([_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_) – The input structure.

  * **indices** (_list_ _[__int_ _] or_ _None_ _,__optional_) – List of indices of atoms allowed to vibrate. Defaults to None.



Returns
    

The vibration calculation result.

Return type
    

[FastVibrationsResult](matlantis_features.features.common.vibration.FastVibrationsResult.html#matlantis_features.features.common.vibration.FastVibrationsResult "matlantis_features.features.common.vibration.FastVibrationsResult")

[ __ previous matlantis_features.features.common.vibration.VibrationFeatureResult ](matlantis_features.features.common.vibration.VibrationFeatureResult.html "previous page") [ next matlantis_features.features.common.vibration.FastVibrationsResult __](matlantis_features.features.common.vibration.FastVibrationsResult.html "next page")
