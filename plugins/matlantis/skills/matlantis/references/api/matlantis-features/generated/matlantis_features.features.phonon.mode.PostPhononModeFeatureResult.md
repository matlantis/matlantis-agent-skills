# matlantis_features.features.phonon.mode.PostPhononModeFeatureResult#

_class _matlantis_features.features.phonon.mode.PostPhononModeFeatureResult(_k_point : Tuple[float, float, float]_, _repeat_of_cell : Tuple[int, int, int]_, _amplitude : float_, _n_images : int_, _mode_indices : List[int]_, _frequencies : List[float]_, _trajectories : List[List[Atoms]]_, _q_point : Optional[Tuple[float, float, float]] = None_)[[source]](../_modules/matlantis_features/features/phonon/mode.html#PostPhononModeFeatureResult)#
    

Bases: `object`

A dataclass for result of PostPhononModeFeature.

Methods

`__init__`(k_point, repeat_of_cell, amplitude, ...) |   
---|---  
`plot`(save_path[, rotation]) | Creates gif files to show the movement of atoms in each phonon mode.  
`to_traj`(save_path) | Create ase.io.Trajectory to save the movement of atoms in each phonon mode.  
  
__init__(_k_point : Tuple[float, float, float]_, _repeat_of_cell : Tuple[int, int, int]_, _amplitude : float_, _n_images : int_, _mode_indices : List[int]_, _frequencies : List[float]_, _trajectories : List[List[Atoms]]_, _q_point : Optional[Tuple[float, float, float]] = None_) → None#
    

plot(_save_path : str_, _rotation : str = '0x,0y,0z'_) → None[[source]](../_modules/matlantis_features/features/phonon/mode.html#PostPhononModeFeatureResult.plot)#
    

Creates gif files to show the movement of atoms in each phonon mode.

Parameters
    

  * **save_path** (_str_) – The directory to which the gif files will be saved. The file will be named as ‘mode’ + band index + ‘.gif’. The band index is given in ascending order of phonon frequency.

  * **rotation** (_str_ _,__optional_) – The viewing angle. Defaults to ‘0x,0y,0z’.




to_traj(_save_path : str_) → None[[source]](../_modules/matlantis_features/features/phonon/mode.html#PostPhononModeFeatureResult.to_traj)#
    

Create ase.io.Trajectory to save the movement of atoms in each phonon mode.

Parameters
    

**save_path** (_str_) – The directory to which the gif files will be saved. The file will be named as ‘mode’ + band index + ‘.traj’. The band index is given in ascending order of phonon frequency.

Attributes

dict#
    

Dictionary representation of the PostPhononModeFeatureResult.

Returns
    

Phonon modes in dictionary format.

Return type
    

dict[str, Any]

q_point _: Optional[Tuple[float, float, float]]__ = None_#
    

k_point _: Tuple[float, float, float]_#
    

repeat_of_cell _: Tuple[int, int, int]_#
    

amplitude _: float_#
    

n_images _: int_#
    

mode_indices _: List[int]_#
    

frequencies _: List[float]_#
    

trajectories _: List[List[Atoms]]_#
    

[ __ previous matlantis_features.features.phonon.mode.PostPhononModeFeature ](matlantis_features.features.phonon.mode.PostPhononModeFeature.html "previous page") [ next matlantis_features.features.phonon.symmetry __](matlantis_features.features.phonon.symmetry.html "next page")
