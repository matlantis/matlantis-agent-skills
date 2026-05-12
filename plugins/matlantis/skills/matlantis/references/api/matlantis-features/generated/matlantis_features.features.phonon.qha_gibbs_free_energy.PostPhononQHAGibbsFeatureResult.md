# matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult#

_class _matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult(_force_constant_list : List[[ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")]_, _unit_cell_atoms : List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _kpts : List[int]_, _temperatures : List[float]_, _pressures : List[float]_, _thermochemistry : List[[PostPhononThermochemistryFeatureResult](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult "matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult")]_, _equilibrium_volume : ndarray_, _helmholtz_free_energy : ndarray_, _gibbs_free_energy : ndarray_)[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#PostPhononQHAGibbsFeatureResult)#
    

Bases: `object`

A dataclass for result of PostPhononQHAGibbsFeatureResult.

Methods

`__init__`(force_constant_list, ...) |   
---|---  
`plot`([plt_name]) | Plots all calculation results from QHA.  
`plot_3d`() | Plots the surface of Gibbs free energy G(T, P).  
  
__init__(_force_constant_list : List[[ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")]_, _unit_cell_atoms : List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _kpts : List[int]_, _temperatures : List[float]_, _pressures : List[float]_, _thermochemistry : List[[PostPhononThermochemistryFeatureResult](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult "matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult")]_, _equilibrium_volume : ndarray_, _helmholtz_free_energy : ndarray_, _gibbs_free_energy : ndarray_) → None#
    

plot(_plt_name : Optional[str] = None_) → Figure[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#PostPhononQHAGibbsFeatureResult.plot)#
    

Plots all calculation results from QHA.

Parameters
    

**plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of Gibbs free energy will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

Returns
    

The plot of Gibbs free energy as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

plot_3d() → Figure[[source]](../_modules/matlantis_features/features/phonon/qha_gibbs_free_energy.html#PostPhononQHAGibbsFeatureResult.plot_3d)#
    

Plots the surface of Gibbs free energy G(T, P).

Returns
    

The plot of Gibbs free energy as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

Attributes

force_constant_list _: List[[ForceConstantFeatureResult](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html#matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult "matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult")]_#
    

unit_cell_atoms _: List[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_#
    

kpts _: List[int]_#
    

temperatures _: List[float]_#
    

pressures _: List[float]_#
    

thermochemistry _: List[[PostPhononThermochemistryFeatureResult](matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult "matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult")]_#
    

equilibrium_volume _: ndarray_#
    

helmholtz_free_energy _: ndarray_#
    

gibbs_free_energy _: ndarray_#
    

[ __ previous matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature ](matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature.html "previous page") [ next matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram __](matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram.html "next page")
