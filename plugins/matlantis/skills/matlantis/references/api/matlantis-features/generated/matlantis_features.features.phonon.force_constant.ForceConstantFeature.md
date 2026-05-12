# matlantis_features.features.phonon.force_constant.ForceConstantFeature#

_class _matlantis_features.features.phonon.force_constant.ForceConstantFeature(_supercell : Tuple[int, int, int]_, _delta : float_, _method : str = 'Frederiksen'_, _direction : str = 'central'_, _symmetrize : int = 3_, _acoustic : bool = True_, _cutoff : Optional[float] = None_, _use_symmetry : bool = False_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/phonon/force_constant.html#ForceConstantFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the force constant for the given structure.

The force constant, which is the basis of phonon analysis, will be calculated with the finite displacement method, in which the second derivative of the potential energy is estimated by displacing each atom with a small distance within an enlarged supercell. The results will be used in the post analysis, such as PostPhononBandFeature.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/force_constant.html).

Methods

`__init__`(supercell, delta[, method, ...]) | Initialize an instance.  
---|---  
`__call__`(atoms[, primitive_matrix, prec]) | Calculates force constant for the input structure.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_supercell : Tuple[int, int, int]_, _delta : float_, _method : str = 'Frederiksen'_, _direction : str = 'central'_, _symmetrize : int = 3_, _acoustic : bool = True_, _cutoff : Optional[float] = None_, _use_symmetry : bool = False_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/phonon/force_constant.html#ForceConstantFeature.__init__)#
    

Initialize an instance.

Parameters
    

  * **supercell** (_tuple_ _[__int_ _,__int_ _,__int_ _]_) – The enlarged supercell that is created by elongating the input structure by (n, m, l) times along the axes.

  * **delta** (_float_) – The distance of finite displacement applied to each atom. The unit is Angstrom.

  * **method** (_str_ _,__optional_) – The method to calculate force constant from the displaced supercells and their forces. This parameter can be ‘Standard’ or ‘Frederiksen’. Defaults to ‘Frederiksen’.

  * **direction** (_str_ _,__optional_) – This parameter discribe how the displacment was found. Must be in ‘forward’, ‘backward’ or ‘central’. When ‘direction=forward’/ ‘direction=backward’, each atom is displaced along the positive/negative direction of the axes. When ‘direction=center’, each atom is displaced along both positive and negative directions of the axes. Defaults to ‘central’.

  * **symmetrize** (_int_ _,__optional_) – After restoring the acoustic sum rule, the symmetry of the force constant is broken. The parameter specifies the number of interactions to symmetrize the force constant. Defaults to 3.

  * **acoustic** (_bool_ _,__optional_) – Restores the acoustic sum rule of the force constant. Defaults to True.

  * **cutoff** (_float_ _or_ _None_ _,__optional_) – The corresponding elements in the force constant are set to 0 if the interatomic distance between atoms is larger than cutoff. If None is provided, the full force constant will be used. Defaults to None.

  * **use_symmetry** (_bool_) – If True, the symmetry of input structure will be used to reduce the number of force estimations. Defaults to False.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.




__call__(_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_, _primitive_matrix : Optional[Union[ndarray, str]] = None_, _prec : float = 0.001_) → [ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")#
    

Calculates force constant for the input structure.

Parameters
    

  * **atoms** ([_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_) – The input structure.

  * **primitive_matrix** (_np.ndarray_ _or_ _str_ _or_ _None_ _,__optional_) – When ‘primitive_matrix’ is not None, the input structure will be transformed to the primitive cell according to the axes provided here. Then, the following force constant calculation will be performed with the reduced primitive cell. Defaults to None.

  * **prec** (_float_ _,__optional_) – Wrapper the input structure into the primitive cell with certain precision. The distance between the wrapped atom and the atom in the primitive cell should be smaller than prec. Defaults to 0.001.



Returns
    

The force constant.

Return type
    

[ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")

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

[ __ previous matlantis_features.features.phonon.force_constant ](matlantis_features.features.phonon.force_constant.html "previous page") [ next matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult __](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html "next page")
