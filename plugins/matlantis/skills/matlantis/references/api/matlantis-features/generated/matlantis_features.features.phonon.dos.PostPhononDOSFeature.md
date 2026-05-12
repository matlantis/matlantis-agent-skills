# matlantis_features.features.phonon.dos.PostPhononDOSFeature#

_class _matlantis_features.features.phonon.dos.PostPhononDOSFeature(_partial : bool = False_)[[source]](../_modules/matlantis_features/features/phonon/dos.html#PostPhononDOSFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the phonon DOS from the force constant.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/phonon_dos.html).

Methods

`__init__`([partial]) | Initialize an instance.  
---|---  
`__call__`(force_constant, kpts[, scheme, ...]) | Calculates the phonon density of states (DOS) from the force constant.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_partial : bool = False_) → None[[source]](../_modules/matlantis_features/features/phonon/dos.html#PostPhononDOSFeature.__init__)#
    

Initialize an instance.

Parameters
    

**partial** (_bool_ _,__optional_) – Whether to calculate the element projected DOS. Defaults to False.

__call__(_force_constant : [ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")_, _kpts : List[int]_, _scheme : str = 'mp'_, _shift : Optional[ndarray] = None_, _freq_min : Optional[float] = None_, _freq_max : Optional[float] = None_, _freq_bin : Optional[float] = None_, _unit : str = 'meV'_) → [PostPhononDOSFeatureResult](matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult "matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult")#
    

Calculates the phonon density of states (DOS) from the force constant.

Parameters
    

  * **force_constant** ([_ForceConstantFeatureResult_](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")) – The calculation result of the ForceConstantFeature.

  * **kpts** (_list_ _[__int_ _]_) – The number of mesh points along each axis.

  * **scheme** (_str_ _,__optional_) – The type of mesh grid: either ‘mp’ or ‘gamma’. The mesh grid is determined by the Monkhorst-Pack scheme (‘mp’) or Gamma-centered scheme (‘gamma’). Defaults to ‘mp’.

  * **shift** (_np.ndarray_ _or_ _None_ _,__optional_) – The shift of the mesh grid in the direction along the corresponding reciprocal axes. Defaults to None.

  * **freq_min** (_float_ _or_ _None_ _,__optional_) – The minimum frequency of phonon DOS calculation range. Defaults to None.

  * **freq_max** (_float_ _or_ _None_ _,__optional_) – The maximum frequency of phonon DOS calculation range. Defaults to None.

  * **freq_bin** (_float_ _or_ _None_ _,__optional_) – The small interval in which the density of states will be calculated. Defaults to None.

  * **unit** (_str_ _,__optional_) – The unit of phonon frequency. Must be one of: ‘meV’, ‘eV’, ‘THz’ or ‘cm^-1’. Defaults to ‘meV’.



Returns
    

The calculation results of phonon DOS.

Return type
    

[PostPhononDOSFeatureResult](matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult "matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult")

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

[ __ previous matlantis_features.features.phonon.dos ](matlantis_features.features.phonon.dos.html "previous page") [ next matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult __](matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html "next page")
