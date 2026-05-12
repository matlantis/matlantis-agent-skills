# matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature#

_class _matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeature[[source]](../_modules/matlantis_features/features/md/post_features/diffusion.html#PostMDDiffusionFeature)#
    

Bases: [`FeatureBase`](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

The matlantis-feature for calculating the diffusion coefficient using the MD trajectory.

Note

The breaking change is made in `PostMDDiffusionFeature` at version 0.8.0. The new implementation is fast, and enables a more comprehensive analysis of diffusion properties. The old version is renamed as [`OldPostMDDiffusionFeature`](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html#matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature "matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature").

References

For more information, please refer to [Matlantis Guidebook](https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/diffusion.html).

Methods

`__init__`() | Initialize an instance.  
---|---  
`__call__`(md_results[, init_time, stride, ...]) | Calculate the diffusion coefficient of each element from the MD trajectory.  
`attach_ctx`([ctx]) | Attach the feature to matlantis_features.utils.Context.  
`check_estimator_fn`(estimator_fn) | Checks if the given estimator function is None and output a warning if so.  
`cost_estimate`([atoms]) | Estimate the cost of the feature.  
`from_dict`(res) | Construct a FeatureBase object from serialized dict.  
`get_savedir_from_ctx`() | Get the temporary save directory from the context.  
`init_scope`() | Context manager that enable to set attribution of the feature.  
`repeat`(n_repeat) | Set the maximum number of times that allowed to run the __call__ function.  
`to_dict`() | Dictionary representation of the FeatureBase.  
  
__init__() → None#
    

Initialize an instance.

__call__(_md_results : [MDFeatureResult](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")_, _init_time : float = 0.0_, _stride : int = 1_, _atom_indices : Optional[List[int]] = None_, _molecule : Union[bool, List[List[int]]] = False_, _number_of_segments : int = 1_, _direction : Optional[ndarray] = None_, _method : str = 'normal'_, _effective_msd_range : Optional[Tuple[float, float]] = None_) → [PostMDDiffusionFeatureResult](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult "matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult")#
    

Calculate the diffusion coefficient of each element from the MD trajectory.

Parameters
    

  * **md_results** ([_MDFeatureResult_](matlantis_features.features.md.md.MDFeatureResult.html#matlantis_features.features.md.md.MDFeatureResult "matlantis_features.features.md.md.MDFeatureResult")) – MD feature’s result object containing the trajectory.

  * **init_time** (_float_ _,__optional_) – Trajectory frames before this time (in fs) is ignored for the calculation. Defaults to 0.0 (i.e., all frames are used for calculation).

  * **stride** (_int_ _,__optional_) – Size of the skip strie of the trajectory frames. Defaults to 1 (i.e., all frames are used for calculation).

  * **atom_indices** (_list_ _[__int_ _] or_ _None_ _,__optional_) – The indices of atoms whose diffusion coefficient is to be calculated explicitly. Defaults to None.

  * **molecule** (_bool_ _or_ _list_ _[__list_ _[__int_ _]__]__,__optional_) – The indices of atoms who form a molecule. If it is used, the diffusion coefficient of the center of mass of molecules will be calculated. The results will be saved in the PostMDDiffusionFeatureResult.diffusion_coefficient_molecule. If True, all the atoms will be regard as one molecule. Defaults to False.

  * **number_of_segments** (_int_ _,__optional_) – Divides the given trajectory in to segments to allow statistical analysis, such as standard deviation of diffuion coefficient. Defaults to 1.

  * **direction** (_np.ndarray_ _or_ _None_ _,__optional_) – Calculate the diffusion coefficient along specific directions. It must be a nx3 array where n is number of directions. The result will be saved in PostMDDiffusionFeatureResult.diffusion_coefficient as, e.g. H_direc_1, O_direc_2 etc. If None, only the diffusion coefficient along x, y and z directions is calculated. Defaults to None

  * **method** (_str_ _,__optional_) – The method to calculate mean square displacement (MSD). Now, two methods, i.e. “normal” and “segment”, are supported. Defaults to “normal”.

  * **effective_msd_range** (_tuple_ _[__float_ _,__float_ _] or_ _None_ _,__optional_) – A tuple of two floats representing the range of MSD sequence to be used in the calculation of diffusion coefficient. The values should be in the range of [0.0, 1.0]. If None, all the mean square displacements will be used. Defaults to None.



Returns
    

Result dataclass object.

Return type
    

[PostMDDiffusionFeatureResult](matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult.html#matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult "matlantis_features.features.md.post_features.diffusion.PostMDDiffusionFeatureResult")

attach_ctx(_ctx : Optional[Context] = None_) → None#
    

Attach the feature to matlantis_features.utils.Context.

Parameters
    

**ctx** (_Context_ _or_ _None_ _,__optional_) – The matlantis_features.utils.Context object. Defaults to None.

check_estimator_fn(_estimator_fn : Optional[Callable[[], Estimator]]_) → None#
    

Checks if the given estimator function is None and output a warning if so.

Parameters
    

**estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail. Defaults to None.

cost_estimate(_atoms : Optional[Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]] = None_) → FeatureCost#
    

Estimate the cost of the feature.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _None_ _,__optional_) – The input atoms. Defaults to None.

Returns
    

The cost of the feature.

Return type
    

FeatureCost

_classmethod _from_dict(_res : Dict[str, Any]_) → [FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")#
    

Construct a FeatureBase object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized FeatureBase from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

[FeatureBase](matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase "matlantis_features.features.base.FeatureBase")

get_savedir_from_ctx() → Path#
    

Get the temporary save directory from the context.

Returns
    

The temporary save directory .

Return type
    

pathlib.Path

init_scope() → Iterator[None]#
    

Context manager that enable to set attribution of the feature.

Returns
    

Init_scope context manager.

Return type
    

Iterator[None]

repeat(_n_repeat : int_) → Self#
    

Set the maximum number of times that allowed to run the __call__ function.

Parameters
    

**n_repeat** (_int_) – The maximum number of repeats.

Returns
    

The feature.

Return type
    

Self

to_dict() → Dict[str, Any]#
    

Dictionary representation of the FeatureBase.

Returns
    

A dict containing a serialized FeatureBase.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature ](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeature.html "previous page") [ next matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult __](matlantis_features.features.md.post_features.diffusion.OldPostMDDiffusionFeatureResult.html "next page")
