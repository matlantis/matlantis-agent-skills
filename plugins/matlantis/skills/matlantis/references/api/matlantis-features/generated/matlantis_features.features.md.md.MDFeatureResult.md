# matlantis_features.features.md.md.MDFeatureResult#

_class _matlantis_features.features.md.md.MDFeatureResult(_traj_path : Optional[str] = None_, _checkpoint_path : Optional[str] = None_, _temp_dir : Optional[TemporaryDirectory] = None_)[[source]](../_modules/matlantis_features/features/md/md.html#MDFeatureResult)#
    

Bases: `object`

A dataclass for result of MDFeature.

Methods

`__init__`([traj_path, checkpoint_path, temp_dir]) |   
---|---  
`from_traj_obj`(traj_file_name) | Get MDFeatureResult from trajectory file.  
`get_traj_obj`(ctxt_save_dir) | Get ASE's trajectory reader object of the MD result.  
`save`([traj_file_name, checkpoint_file_name]) | Save the MD trajectory and checkpoint files to a different location.  
  
__init__(_traj_path : Optional[str] = None_, _checkpoint_path : Optional[str] = None_, _temp_dir : Optional[TemporaryDirectory] = None_) → None#
    

_classmethod _from_traj_obj(_traj_file_name : str_) → MDFeatureResult[[source]](../_modules/matlantis_features/features/md/md.html#MDFeatureResult.from_traj_obj)#
    

Get MDFeatureResult from trajectory file.

Parameters
    

**traj_file_name** (_str_) – Path to the trajectory file.

Returns
    

MDFeatureResult.

Return type
    

MDFeatureResult

get_traj_obj(_ctxt_save_dir : Path_) → TrajectoryReader[[source]](../_modules/matlantis_features/features/md/md.html#MDFeatureResult.get_traj_obj)#
    

Get ASE’s trajectory reader object of the MD result.

Parameters
    

**ctxt_save_dir** (_Path_) – Context’s temporary directory location. Ignored if this obj contains trajectory path in absolute path.

Returns
    

TrajectoryReader object for this MD result.

Return type
    

TrajectoryReader

save(_traj_file_name : Optional[str] = None_, _checkpoint_file_name : Optional[str] = None_) → None[[source]](../_modules/matlantis_features/features/md/md.html#MDFeatureResult.save)#
    

Save the MD trajectory and checkpoint files to a different location.

Parameters
    

  * **traj_file_name** (_str_ _or_ _None_ _,__optional_) – The new location of MD trajectory file. Defaults to None.

  * **checkpoint_file_name** (_str_ _or_ _None_ _,__optional_) – The new location of MD checkpoint file. Defaults to None.




Attributes

checkpoint_path _: Optional[str]__ = None_#
    

temp_dir _: Optional[TemporaryDirectory]__ = None_#
    

traj_path _: Optional[str]__ = None_#
    

[ __ previous matlantis_features.features.md.md.MDFeature ](matlantis_features.features.md.md.MDFeature.html "previous page") [ next matlantis_features.features.md.md_extensions __](matlantis_features.features.md.md_extensions.html "next page")
