# matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult#

_class _matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult(_temperature : Dict[str, List[float]]_, _pressure : Dict[str, List[float]]_, _formation_enthalpy : Dict[str, List[float]]_)[[source]](../_modules/matlantis_features/features/common/gas_formation_enthalpy.html#ComplexGasFormationEnthalpyFeatureResult)#
    

Bases: `object`

A dataclass for result of ComplexFormationEnthalpyFeature.

Methods

`__init__`(temperature, pressure, ...) |   
---|---  
`to_dataframe`([csv_name]) | Converts the formation enthalpy calculation result to pandas.DataFrame.  
  
__init__(_temperature : Dict[str, List[float]]_, _pressure : Dict[str, List[float]]_, _formation_enthalpy : Dict[str, List[float]]_) → None#
    

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/common/gas_formation_enthalpy.html#ComplexGasFormationEnthalpyFeatureResult.to_dataframe)#
    

Converts the formation enthalpy calculation result to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The csv file path to which the dataframe will be saved. If ‘None’ is provided, the dataframe will not be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains temperature, pressure and formation_enthalpy.

Return type
    

pd.DataFrame

Attributes

dict#
    

Dictionary representation of the ComplexFormationEnthalpyFeatureResult.

Returns
    

Formation enthalpy in dictionary format.

Return type
    

dict[str, Any]

temperature _: Dict[str, List[float]]_#
    

pressure _: Dict[str, List[float]]_#
    

formation_enthalpy _: Dict[str, List[float]]_#
    

[ __ previous matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature ](matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature.html "previous page") [ next matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeatureResult __](matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeatureResult.html "next page")
