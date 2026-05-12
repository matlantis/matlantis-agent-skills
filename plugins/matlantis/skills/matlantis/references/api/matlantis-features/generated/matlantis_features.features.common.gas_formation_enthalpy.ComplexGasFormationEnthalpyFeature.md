# matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature#

_class _matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature(_opt_feature : [OptFeatureBase](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")_, _vibration_feature : [VibrationFeature](matlantis_features.features.common.vibration.VibrationFeature.html#matlantis_features.features.common.vibration.VibrationFeature "matlantis_features.features.common.vibration.VibrationFeature")_)[[source]](../_modules/matlantis_features/features/common/gas_formation_enthalpy.html#ComplexGasFormationEnthalpyFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating gas formation enthalpy.

Methods

`__init__`(opt_feature, vibration_feature) | Initialize an instance.  
---|---  
`__call__`(atoms, temperatures, pressures) | Calculate the formation enthalpy of the gas molecule.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__(_opt_feature : [OptFeatureBase](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")_, _vibration_feature : [VibrationFeature](matlantis_features.features.common.vibration.VibrationFeature.html#matlantis_features.features.common.vibration.VibrationFeature "matlantis_features.features.common.vibration.VibrationFeature")_) → None[[source]](../_modules/matlantis_features/features/common/gas_formation_enthalpy.html#ComplexGasFormationEnthalpyFeature.__init__)#
    

Initialize an instance.

Note

Since `ComplexGasFormationEnthalpyFeature` instantiates an estimator in __init__, it won’t work with multiprocessing parallelization. If you want to use this feature in multiple processes, please make sure you use loky for the backend. For more information, please refer to the [Matlantis guidebook](/api/resource/documents/matlantis-guidebook/en/notebooks.html#the-use-of-estimator-in-multi-process-computations).

Parameters
    

  * **opt_feature** ([_OptFeatureBase_](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")) – Optimizer feature object used for structure optimization before the enthalpy calculation.

  * **vibration_feature** ([_VibrationFeature_](matlantis_features.features.common.vibration.VibrationFeature.html#matlantis_features.features.common.vibration.VibrationFeature "matlantis_features.features.common.vibration.VibrationFeature")) – Feature object for calculating molecule vibration.




__call__(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _temperatures : List[float]_, _pressures : List[float]_) → [ComplexGasFormationEnthalpyFeatureResult](matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult "matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult")#
    

Calculate the formation enthalpy of the gas molecule.

Parameters
    

  * **atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – The input structure.

  * **temperatures** (_list_ _[__float_ _]_) – The temperatures at which the formation enthalpy will be calculated.

  * **pressures** (_list_ _[__float_ _]_) – The pressures at which the formation enthalpy will be calculated. The unit of pressure is eV/Angstrom^3.



Returns
    

The calculation results of formation enthalpy.

Return type
    

[ComplexGasFormationEnthalpyFeatureResult](matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult "matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult")

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

[ __ previous matlantis_features.features.common.gas_formation_enthalpy ](matlantis_features.features.common.gas_formation_enthalpy.html "previous page") [ next matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature __](matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature.html "next page")
