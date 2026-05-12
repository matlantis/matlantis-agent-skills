# matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult#

_class _matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult(_temperature : Dict[str, List[float]]_, _pressure : Dict[str, List[float]]_, _zero_point_energy : float_, _enthalpy : Dict[str, List[float]]_, _entropy : Dict[str, List[float]]_, _gibbs_free_energy : Dict[str, List[float]]_)[[source]](../_modules/matlantis_features/features/common/gas_thermo.html#PostVibrationGasThermoFeatureResult)#
    

Bases: `object`

A dataclass for result of PostVibrationGasThermoFeature.

Methods

`__init__`(temperature, pressure, ...) |   
---|---  
`to_dataframe`([csv_name]) | Converts the thermochemical properties to pandas.DataFrame.  
  
__init__(_temperature : Dict[str, List[float]]_, _pressure : Dict[str, List[float]]_, _zero_point_energy : float_, _enthalpy : Dict[str, List[float]]_, _entropy : Dict[str, List[float]]_, _gibbs_free_energy : Dict[str, List[float]]_) → None#
    

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/common/gas_thermo.html#PostVibrationGasThermoFeatureResult.to_dataframe)#
    

Converts the thermochemical properties to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The csv file path to which the dataframe will be saved. If ‘None’ is provided, the dataframe will not be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains temperature, pressure, entropy, enthalpy and helmholtz_free_energy,

Return type
    

pd.DataFrame

Attributes

dict#
    

Dictionary representation of the PostVibrationGasThermoFeatureResult.

Returns
    

Thermochemistry properties in dictionary format.

Return type
    

dict[str, Any]

temperature _: Dict[str, List[float]]_#
    

pressure _: Dict[str, List[float]]_#
    

zero_point_energy _: float_#
    

enthalpy _: Dict[str, List[float]]_#
    

entropy _: Dict[str, List[float]]_#
    

gibbs_free_energy _: Dict[str, List[float]]_#
    

[ __ previous matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeature ](matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeature.html "previous page") [ next Gas Formation Enthalpy Features __](../matlantis_features.features.common.gas_formation_enthalpy.html "next page")
