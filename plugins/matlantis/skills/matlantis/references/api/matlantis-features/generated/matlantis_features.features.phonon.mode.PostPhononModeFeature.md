# matlantis_features.features.phonon.mode.PostPhononModeFeature#

_class _matlantis_features.features.phonon.mode.PostPhononModeFeature[[source]](../_modules/matlantis_features/features/phonon/mode.html#PostPhononModeFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature to visualize how atoms move under a certain phonon mode.

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(force_constant[, k_point, ...]) | Calculates the phonon mode from the force constant.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None[[source]](../_modules/matlantis_features/features/phonon/mode.html#PostPhononModeFeature.__init__)#
    

Initialize an instance.

__call__(_force_constant : [ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")_, _k_point : Tuple[float, float, float] = (0.0, 0.0, 0.0)_, _repeat_of_cell : Tuple[int, int, int] = (3, 3, 3)_, _mode_indices : Optional[Union[int, List[int]]] = None_, _n_images : int = 20_, _amplitude : float = 5.0_, _q_point : Optional[Tuple[float, float, float]] = None_) → [PostPhononModeFeatureResult](matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html#matlantis_features.features.phonon.mode.PostPhononModeFeatureResult "matlantis_features.features.phonon.mode.PostPhononModeFeatureResult")#
    

Calculates the phonon mode from the force constant.

Parameters
    

  * **force_constant** ([_ForceConstantFeatureResult_](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")) – The force contant result obtained from matlantis_features.features.phonon.RunForceConstantFeature

  * **k_point** (_tuple_ _[__float_ _,__float_ _,__float_ _]__,__optional_) – The k-point in which the phonon mode will be calculated and the atom movement will be viewed. Defaults to [0., 0., 0.].

  * **repeat_of_cell** (_tuple_ _[__int_ _,__int_ _,__int_ _]__,__optional_) – To view the phonon mode at a large supercell. Defaults to [3, 3, 3].

  * **mode_indices** (_int_ _or_ _list_ _[__int_ _] or_ _None_ _,__optional_) – The index of branches. If None is provided, the phonon modes will be calculated in all branches.

  * **n_images** (_int_ _,__optional_) – Number of images in one period. Defaults to 20.

  * **amplitude** (_float_ _,__optional_) – The amplitude of vibration. Defaults to 1.0.

  * **q_point** (_tuple_ _[__float_ _,__float_ _,__float_ _] or_ _None_ _,__optional_) – Deprecated parameter. The q-vector in which the phonon mode will be calculated and the atom movement will be viewed. Defaults to [0., 0., 0.].



Returns
    

A set of structure images that show the periodic movement of each atom under phonon modes of given a k-point.

Return type
    

[PostPhononModeFeatureResult](matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html#matlantis_features.features.phonon.mode.PostPhononModeFeatureResult "matlantis_features.features.phonon.mode.PostPhononModeFeatureResult")

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

[ __ previous matlantis_features.features.phonon.mode ](matlantis_features.features.phonon.mode.html "previous page") [ next matlantis_features.features.phonon.mode.PostPhononModeFeatureResult __](matlantis_features.features.phonon.mode.PostPhononModeFeatureResult.html "next page")
