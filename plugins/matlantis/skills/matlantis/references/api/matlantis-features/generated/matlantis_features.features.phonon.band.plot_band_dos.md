# matlantis_features.features.phonon.band.plot_band_dos#

matlantis_features.features.phonon.band.plot_band_dos(_band_results : [PostPhononBandFeatureResult](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult "matlantis_features.features.phonon.band.PostPhononBandFeatureResult")_, _dos_results : [PostPhononDOSFeatureResult](matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult "matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult")_, _plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/phonon/band.html#plot_band_dos)#
    

Plot the phonon band structure and phonon density of states in one figure.

Parameters
    

  * **band_results** ([_PostPhononBandFeatureResult_](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult "matlantis_features.features.phonon.band.PostPhononBandFeatureResult")) – The result of PostPhononBandFeature.

  * **dos_results** ([_PostPhononDOSFeatureResult_](matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult "matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult")) – The result of PostPhononDOSFeature.

  * **plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the figure will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.



Returns
    

The plot of the phonon band structure and phonon density of states as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

[ __ previous matlantis_features.features.phonon.band.PostPhononBandFeatureResult ](matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html "previous page") [ next matlantis_features.features.phonon.dos __](matlantis_features.features.phonon.dos.html "next page")
