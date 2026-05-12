# matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult#

_class _matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult(_diffusion_coefficient : Dict[str, Dict[str, float]]_, _diffusion_coefficient_std : Dict[str, Dict[str, float]]_, _diffusion_coefficient_molecule : Dict[str, Dict[str, float]]_, _diffusion_coefficient_molecule_std : Dict[str, Dict[str, float]]_, _mean_squared_displacement : Dict[str, ndarray]_, _mean_squared_displacement_molecule : Dict[str, ndarray]_, _timestep : float_, _init_time : float_, _stride : int_, _atom_indices : Optional[List[int]]_, _molecule : Union[bool, List[List[int]]]_, _number_of_segments : int_, _direction : Optional[ndarray]_, _method : str_, _effective_msd_range : Optional[Tuple[float, float]]_)[[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#PostMDDiffusionFeatureResult)#
    

Bases: `object`

A dataclass for result of PostMDDiffusionFeature.

Methods

`__init__`(diffusion_coefficient, ...) |   
---|---  
`plot`([plt_name, molecule, show_fit_line, ...]) | Plot the mean square displacement.  
  
__init__(_diffusion_coefficient : Dict[str, Dict[str, float]]_, _diffusion_coefficient_std : Dict[str, Dict[str, float]]_, _diffusion_coefficient_molecule : Dict[str, Dict[str, float]]_, _diffusion_coefficient_molecule_std : Dict[str, Dict[str, float]]_, _mean_squared_displacement : Dict[str, ndarray]_, _mean_squared_displacement_molecule : Dict[str, ndarray]_, _timestep : float_, _init_time : float_, _stride : int_, _atom_indices : Optional[List[int]]_, _molecule : Union[bool, List[List[int]]]_, _number_of_segments : int_, _direction : Optional[ndarray]_, _method : str_, _effective_msd_range : Optional[Tuple[float, float]]_) → None#
    

plot(_plt_name : Optional[str] = None_, _molecule : bool = False_, _show_fit_line : bool = False_, _show_effective_range : bool = False_) → Figure[[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#PostMDDiffusionFeatureResult.plot)#
    

Plot the mean square displacement.

Parameters
    

  * **plt_name** (_str_ _or_ _None_ _,__optional_) – File name to write the plot. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

  * **molecule** (_bool_ _,__optional_) – Plot the mean sqare displacement of molecular center of mass. Defaults to None.

  * **show_fit_line** (_bool_ _,__optional_) – Plot the fitting line. Defaults to False.

  * **show_effective_range** (_bool_ _,__optional_) – Plot the effective range of msd that used for the linear fitting. Defaults to False.



Returns
    

Resulting plotly’s graph object.

Return type
    

go.Figure

Attributes

dict#
    

Convert the result to the dict object.

Returns
    

Resulting dict object.

Return type
    

dict[str, Any]

diffusion_coefficient _: Dict[str, Dict[str, float]]_#
    

diffusion_coefficient_std _: Dict[str, Dict[str, float]]_#
    

diffusion_coefficient_molecule _: Dict[str, Dict[str, float]]_#
    

diffusion_coefficient_molecule_std _: Dict[str, Dict[str, float]]_#
    

mean_squared_displacement _: Dict[str, ndarray]_#
    

mean_squared_displacement_molecule _: Dict[str, ndarray]_#
    

timestep _: float_#
    

init_time _: float_#
    

stride _: int_#
    

atom_indices _: Optional[List[int]]_#
    

molecule _: Union[bool, List[List[int]]]_#
    

number_of_segments _: int_#
    

direction _: Optional[ndarray]_#
    

method _: str_#
    

effective_msd_range _: Optional[Tuple[float, float]]_#
    

[ __ previous matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult ](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.html "previous page") [ next matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients __](matlantis_features.features.md.post_features.diffusion.fit_diffusion_coefficients.html "next page")
