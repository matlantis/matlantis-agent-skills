# matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature#

_class _matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature[[source]](../_modules/matlantis_features/features/phonon/thermo.html#PostPhononThermochemistryFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating thermochemistry properties from phonon.

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(force_constant, kpts[, tmin, tmax, ...]) | Calculates the thermochemistry properties from the force constant.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None[[source]](../_modules/matlantis_features/features/phonon/thermo.html#PostPhononThermochemistryFeature.__init__)#
    

Initialize an instance.

__call__(_force_constant : [ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")_, _kpts : List[int]_, _tmin : float = 0.0_, _tmax : float = 1000.0_, _tstep : float = 10.0_, _temperatures : Optional[List[float]] = None_, _include_lattice_energy : bool = False_, _scheme : str = 'mp'_, _shift : Optional[ndarray] = None_, _cutoff_frequency : float = 0.0_) → [PostPhononThermochemistryFeatureResult](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult "matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult")#
    

Calculates the thermochemistry properties from the force constant.

Parameters
    

  * **force_constant** ([_ForceConstantFeatureResult_](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")) – The calculation result of the ForceConstantFeature.

  * **kpts** (_list_ _[__int_ _]_) – The number of mesh points along each axis.

  * **tmin** (_float_ _,__optional_) – The lower limit of temperate range to be calculated. Defaults to 0.0.

  * **tmax** (_float_ _,__optional_) – The upper limit of temperate range to be calculated. Defaults to 1000.0.

  * **tstep** (_float_ _,__optional_) – Creates a evenly spaced set of temperatures in the range of [tmin, tmax) at intervals of tstep, at which the thermochemical properties will be calculated. Defaults to 10.0.

  * **temperatures** (_list_ _[__float_ _] or_ _None_ _,__optional_) – The temperatures at which the thermochemistry properties will be calculated. If None is provided, the calculation temperatures will be generated according to ‘tmin’, ‘tmax’ and ‘tstep’. Defaults to None.

  * **include_lattice_energy** (_bool_ _,__optional_) – Add the potential energy of the input structure to the result of internal energy and helmholtz free energy. Defaults to False.

  * **scheme** (_str_ _,__optional_) – The type of mesh grid: either ‘mp’ or ‘gamma’. The mesh grid is determined by the Monkhorst-Pack scheme (‘mp’) or Gamma-centered scheme (‘gamma’). Defaults to ‘mp’.

  * **shift** (_np.ndarray_ _or_ _None_ _,__optional_) – An array of shape (3,) containing the shift of the mesh grid in the direction along the corresponding reciprocal axes. Defaults to None.

  * **cutoff_frequency** (_float_ _,__optional_) – The vibration modes with frequency smaller than ‘cutoff_frequency’ will not be used in the thermochemical property calculation. Imaginary frequencies will be regarded as negative. Defaults to 0.0.



Returns
    

The calculation results of thermochemical properties from phonons.

Return type
    

[PostPhononThermochemistryFeatureResult](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult "matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult")

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

[ __ previous matlantis_features.features.phonon.thermo ](matlantis_features.features.phonon.thermo.html "previous page") [ next matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult __](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html "next page")
