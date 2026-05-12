# matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram#

matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram(_labels : List[str]_, _gibbs_free_energies : List[[PostPhononQHAGibbsFeatureResult](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult "matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult")]_, _interval : int = 501_) → Figure[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#draw_phase_diagram)#
    

Draw the phase diagram from the Gibbs free energy obtained from PostPhononQHAGibbsFeature.

Parameters
    

  * **labels** (_list_ _[__str_ _]_) – The name of each phase.

  * **gibbs_free_energies** (_List_ _[_[_PostPhononQHAGibbsFeatureResult_](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult "matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult") _]_) – The Gibbs free energy calculation results for each phase.

  * **interval** (_int_ _,__optional_) – Number of inteploate points in both temperature and pressure. Defaults to 501.



Returns
    

The phase diagram.

Return type
    

plt.Figure

[ __ previous matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult ](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html "previous page") [ next matlantis_features.features.phonon.utils __](matlantis_features.features.phonon.utils.html "next page")
