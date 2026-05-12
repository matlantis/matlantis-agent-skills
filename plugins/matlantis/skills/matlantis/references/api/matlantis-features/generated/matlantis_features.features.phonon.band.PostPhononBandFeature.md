# matlantis_features.features.phonon.band.PostPhononBandFeature#

_class _matlantis_features.features.phonon.band.PostPhononBandFeature[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the phonon band structure.

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/phonon_band.html).

Methods

`__init__`() | Initiate an instance.  
---|---  
`__call__`(force_constant[, labels, ...]) | Calculates the phonon band structure from the force constant.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeature.__init__)#
    

Initiate an instance.

__call__(_force_constant : [ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")_, _labels : Optional[List[str]] = None_, _special_kpts : Optional[ndarray] = None_, _n_kpts : Optional[List[int]] = None_, _total_n_kpts : int = 100_, _unit : str = 'meV'_, _style : str = 'ase'_) → [PostPhononBandFeatureResult](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult "matlantis_features.features.phonon.band.PostPhononBandFeatureResult")#
    

Calculates the phonon band structure from the force constant.

Parameters
    

  * **force_constant** ([_ForceConstantFeatureResult_](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")) – The calculation result of the ForceConstantFeature.

  * **labels** (_list_ _[__str_ _] or_ _None_ _,__optional_) – The labels of special points that form the band path. Please use ‘,’ ‘|’ or ‘/’ to represent a discontinuous jump. If None is provided, the band path will be automatically generated according to the type of bravais lattice. Defaults to None.

  * **special_kpts** (_np.ndarray_ _or_ _None_ _,__optional_) – The coordinates of special points in reciprocal space. The ‘special_kpts’ should be a numpy.ndarray in the shape of (N, 3), where N is the same as the length of ‘labels’. If ‘None’ is provided, the automatically generated special points will be used. Defaults to None.

  * **n_kpts** (_list_ _[__int_ _] or_ _None_ _,__optional_) – The number of interpolated points in each segment of the path. If None is provided, the values will be estimated from ‘total_n_kpts’. Defaults to None.

  * **total_n_kpts** (_int_ _,__optional_) – The number of interpolated points along the whole band path. This parameter will be ignored if ‘n_kpts’ is specified. Defaults to 100.

  * **unit** (_str_ _,__optional_) – The unit of phonon frequency. Must be one of: ‘meV’, ‘eV’, ‘THz’ or ‘cm^-1’. Defaults to ‘meV’.

  * **style** (_str_ _,__optional_) – If ‘labels’ and ‘special_kpts’ are all not provided, the default k-point path will be used. The parameter ‘style’ can control which kind of default path to use. Currently, ‘ase’, which is ASE style, and ‘phonon_db’, which follows the Phonon database (<http://phonondb.mtl.kyoto-u.ac.jp/>), are supported. Defaults to “ase”.



Returns
    

The calculation results of the phonon band structure.

Return type
    

[PostPhononBandFeatureResult](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult "matlantis_features.features.phonon.band.PostPhononBandFeatureResult")

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

[ __ previous matlantis_features.features.phonon.band ](matlantis_features.features.phonon.band.html "previous page") [ next matlantis_features.features.phonon.band.PostPhononBandFeatureResult __](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html "next page")
