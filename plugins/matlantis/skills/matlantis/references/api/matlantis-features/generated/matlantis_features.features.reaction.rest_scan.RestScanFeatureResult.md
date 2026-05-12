# matlantis_features.features.reaction.rest_scan.RestScanFeatureResult#

_class _matlantis_features.features.reaction.rest_scan.RestScanFeatureResult(_trajectory_logger : TrajectoryLogger_, _traj_path : Optional[str]_, _converged : bool_, _warnflag : Optional[int]_)[[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeatureResult)#
    

Bases: `object`

A dataclass for result of RestScanFeature.

Methods

`__init__`(trajectory_logger, traj_path, ...) |   
---|---  
`from_trajectory`(traj_file_name) | Get RestScanFeatureResult from trajectory file.  
`get_images`([kink_dx, local_minima_de]) | Summerize the images of the restscan result.  
`plot`([plt_name, kink_dx, local_minima_de]) | Plot the energy profile along the images of the restscan result.  
  
__init__(_trajectory_logger : TrajectoryLogger_, _traj_path : Optional[str]_, _converged : bool_, _warnflag : Optional[int]_) → None#
    

_classmethod _from_trajectory(_traj_file_name : Union[str, Path]_) → RestScanFeatureResult[[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeatureResult.from_trajectory)#
    

Get RestScanFeatureResult from trajectory file.

Parameters
    

**traj_file_name** (_str_ _or_ _pathlib.Path_) – Path to the trajectory file.

Returns
    

RestScanFeatureResult.

Return type
    

RestScanFeatureResult

get_images(_kink_dx : float = 0.1_, _local_minima_de : Optional[float] = None_) → List[Atoms][[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeatureResult.get_images)#
    

Summerize the images of the restscan result.

Parameters
    

  * **kink_dx** (_float_ _,__optional_) – Images whose maximum distance is not more than dx are regarded as identical and omitted. Defaults to 0.1.

  * **local_minima_de** (_float_ _or_ _None_ _,__optional_) – If de is given, returns a summary up to the last local minimum. Here, de is the threshold to be regarded as local minimum. Defaults to None.



Returns
    

A list of images.

Return type
    

List[Atoms]

plot(_plt_name : Optional[str] = None_, _kink_dx : float = 0.1_, _local_minima_de : Optional[float] = None_) → Figure[[source]](../_modules/matlantis_features/features/reaction/rest_scan.html#RestScanFeatureResult.plot)#
    

Plot the energy profile along the images of the restscan result.

Parameters
    

  * **plt_name** (_str_ _or_ _None_ _,__optional_) – The file name to which the plot of the energy profile will be saved. Supported formats include ‘png’ ‘jpg’ ‘jpeg’ ‘webp’ ‘svg’ and ‘pdf’. If ‘None’ is provided, the figure will not be saved. Defaults to None.

  * **kink_dx** (_float_ _,__optional_) – Images whose maximum distance is not more than dx are regarded as identical and omitted. Defaults to 0.1.

  * **local_minima_de** (_float_ _or_ _None_ _,__optional_) – If de is given, returns a summary up to the last local minimum. Here, de is the threshold to be regarded as local minimum. Defaults to None.



Returns
    

The plot of the energy profile in the restscan calculation result as a plotly.graph_objects.Figure instance.

Return type
    

go.Figure

Attributes

trajectory_logger _: TrajectoryLogger_#
    

traj_path _: Optional[str]_#
    

converged _: bool_#
    

warnflag _: Optional[int]_#
    

[ __ previous matlantis_features.features.reaction.rest_scan ](matlantis_features.features.reaction.rest_scan.html "previous page") [ next matlantis_features.atoms __](../matlantis_features.atoms.html "next page")
