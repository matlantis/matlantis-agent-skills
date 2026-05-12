# matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult#

_class _matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult(_energy : Dict[str, List[float]]_, _dos : Dict[str, List[float]]_, _unit : str_)[[source]](../_modules/matlantis_features/features/phonon/dos.html#PostPhononDOSFeatureResult)#
    

Bases: `object`

A dataclass for result of PostPhononDOSFeature.

Methods

`__init__`(energy, dos, unit) |   
---|---  
`plot`([plt_name, reverse_axis]) | Plots the phonon density of states (DOS).  
`to_dataframe`([csv_name]) | Convert the phonon DOS to pandas.DataFrame.  
  
__init__(_energy : Dict[str, List[float]]_, _dos : Dict[str, List[float]]_, _unit : str_) → None#
    

plot(_plt_name : Optional[str] = None_, _reverse_axis : bool = False_) → Figure[[source]](../_modules/matlantis_features/features/phonon/dos.html#PostPhononDOSFeatureResult.plot)#
    

Plots the phonon density of states (DOS).

Parameters
    

  * **plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of phonon DOS will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

  * **reverse_axis** (_bool_ _,__optional_) – 

If ‘False’, the x-axis will be ‘energy’ and the y-axis will be density of states. If ‘True’, the x-axis and y-axis will be exchanged.

Defaults to False.



Returns
    

The plot of phonon DOS as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/phonon/dos.html#PostPhononDOSFeatureResult.to_dataframe)#
    

Convert the phonon DOS to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The csv file path to which the dataframe will be saved. If ‘None’ is provided, the dataframe will not be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains the phonon vibration energy and corresponding density of states. The element projected DOS will also be saved if the system contains more than one element.

Return type
    

pd.DataFrame

Attributes

dict#
    

Dictionary representation of the PostPhononDOSFeatureResult.

Returns
    

Phonon DOS in dictionary format.

Return type
    

dict[str, Any]

energy _: Dict[str, List[float]]_#
    

dos _: Dict[str, List[float]]_#
    

unit _: str_#
    

[ __ previous matlantis_features.features.phonon.dos.PostPhononDOSFeature ](matlantis_features.features.phonon.dos.PostPhononDOSFeature.html "previous page") [ next matlantis_features.features.phonon.force_constant __](matlantis_features.features.phonon.force_constant.html "next page")
