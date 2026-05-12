# matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult#

_class _matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult(_kpts : List[int]_, _num_modes : int_, _num_effective_modes : int_, _include_lattice_energy : bool_, _temperature : Dict[str, List[float]]_, _entropy : Dict[str, List[float]]_, _specific_heat : Dict[str, List[float]]_, _internal_energy : Dict[str, List[float]]_, _helmholtz_free_energy : Dict[str, List[float]]_)[[source]](../_modules/matlantis_features/features/phonon/thermo.html#PostPhononThermochemistryFeatureResult)#
    

Bases: `object`

A dataclass for result of PostPhononThermochemistryFeature.

Methods

`__init__`(kpts, num_modes, ...) |   
---|---  
`plot`([plt_name]) | Plots the thermochemistry calculation result from phonon.  
`to_dataframe`([csv_name]) | Converts the thermochemical properties to pandas.DataFrame.  
  
__init__(_kpts : List[int]_, _num_modes : int_, _num_effective_modes : int_, _include_lattice_energy : bool_, _temperature : Dict[str, List[float]]_, _entropy : Dict[str, List[float]]_, _specific_heat : Dict[str, List[float]]_, _internal_energy : Dict[str, List[float]]_, _helmholtz_free_energy : Dict[str, List[float]]_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/phonon/thermo.html#PostPhononThermochemistryFeatureResult.plot)#
    

Plots the thermochemistry calculation result from phonon.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of thermochemistry properties will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

The plot of thermochemical properties as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/phonon/thermo.html#PostPhononThermochemistryFeatureResult.to_dataframe)#
    

Converts the thermochemical properties to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The csv file path to which the dataframe will be saved. If ‘None’ is provided, the dataframe will not be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains temperature, specific_heat, entropy, internal_energy and helmholtz_free_energy,

Return type
    

pd.DataFrame

Attributes

dict#
    

Dictionary representation of the PostPhononThermochemistryFeatureResult.

Returns
    

Thermochemistry properties in dictionary format.

Return type
    

dict[str, Any]

kpts _: List[int]_#
    

num_modes _: int_#
    

num_effective_modes _: int_#
    

include_lattice_energy _: bool_#
    

temperature _: Dict[str, List[float]]_#
    

entropy _: Dict[str, List[float]]_#
    

specific_heat _: Dict[str, List[float]]_#
    

internal_energy _: Dict[str, List[float]]_#
    

helmholtz_free_energy _: Dict[str, List[float]]_#
    

[ __ previous matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature ](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature.html "previous page") [ next matlantis_features.features.phonon.qha_gibbs_free_energy __](matlantis_features.features.phonon.qha_gibbs_free_energy.html "next page")
