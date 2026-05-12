# matlantis_features.features.phonon.band.PostPhononBandFeatureResult#

_class _matlantis_features.features.phonon.band.PostPhononBandFeatureResult(_kpts : ndarray_, _coords : ndarray_, _special_kpt_indices : ndarray_, _labels : List[str]_, _frequency : ndarray_, _unit : str_, _unit_cell_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_)[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeatureResult)#
    

Bases: `object`

A dataclass for result of PostPhononBandFeature.

Methods

`__init__`(kpts, coords, special_kpt_indices, ...) |   
---|---  
`plot`([plt_name]) | Plots the phonon band structure.  
`to_dataframe`([csv_name]) | Converts the phonon band structure to pandas.DataFrame.  
`visual_k_path`([show_reciprocal_basis, ...]) | To visualize the k point path in the 3D plot of the first brillouin zone.  
  
__init__(_kpts : ndarray_, _coords : ndarray_, _special_kpt_indices : ndarray_, _labels : List[str]_, _frequency : ndarray_, _unit : str_, _unit_cell_atoms : [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeatureResult.plot)#
    

Plots the phonon band structure.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of the phonon band structure will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

The plot of the phonon band structure as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

to_dataframe(_csv_name : Optional[str] = None_) → DataFrame[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeatureResult.to_dataframe)#
    

Converts the phonon band structure to pandas.DataFrame.

Parameters
    

**csv_name** (_str_ _or_ _None_ _,__optional_) – The csv file path to which the dataframe will be saved. If ‘None’ is provided, the dataframe will not be saved. Defaults to None.

Returns
    

The pandas.DataFrame which contains the coordinates of the band path, the labels of special points and the phonon vibration frequencies.

Return type
    

pd.DataFrame

visual_k_path(_show_reciprocal_basis : bool = False_, _show_special_k_points : bool = True_, _marker_size : int = 2_) → Figure[[source]](../_modules/matlantis_features/features/phonon/band.html#PostPhononBandFeatureResult.visual_k_path)#
    

To visualize the k point path in the 3D plot of the first brillouin zone.

Parameters
    

  * **show_reciprocal_basis** (_bool_ _,__optional_) – Show the basis vector of the reciprocal cell or not. Defaults to False.

  * **show_special_k_points** (_bool_ _,__optional_) – Show the special k points or not. Defaults to True.

  * **marker_size** (_int_ _,__optional_) – The marker size of the k points. Defaults to 2.



Returns
    

The 3D plot of the k point path.

Return type
    

go.Figure

Attributes

dict#
    

Dictionary representation of the PostPhononBandFeatureResult.

Returns
    

Phonon band structure in dictionary format.

Return type
    

dict[str, Any]

kpts _: ndarray_#
    

coords _: ndarray_#
    

special_kpt_indices _: ndarray_#
    

labels _: List[str]_#
    

frequency _: ndarray_#
    

unit _: str_#
    

unit_cell_atoms _: [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")_#
    

[ __ previous matlantis_features.features.phonon.band.PostPhononBandFeature ](matlantis_features.features.phonon.band.PostPhononBandFeature.html "previous page") [ next matlantis_features.features.phonon.band.plot_band_dos __](matlantis_features.features.phonon.band.plot_band_dos.html "next page")
