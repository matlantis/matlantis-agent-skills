# matlantis_features.features.elasticity.elastic_tensor.stress_from_deform#

matlantis_features.features.elasticity.elastic_tensor.stress_from_deform(_dmat : [DeformTensor](matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor "matlantis_features.features.elasticity.tensors.DeformTensor")_, _atoms : Atoms_, _opt_fn : Callable[[Atoms], Atoms]_, _show_logger : bool_, _pbar : Optional[tqdm_asyncio]_, _estimator_fn : Optional[Callable[[], Estimator]] = None_) → [StressTensor](matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor "matlantis_features.features.elasticity.tensors.StressTensor")[[source]](../_modules/matlantis_features/features/elasticity/elastic_tensor.html#stress_from_deform)#
    

Calculate the stress tensor from the target atoms and deformation tensor, using specified optimization method.

Parameters
    

  * **dmat** ([_DeformTensor_](matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor "matlantis_features.features.elasticity.tensors.DeformTensor")) – Deformation tensor applied to the system.

  * **atoms** (_Atoms_) – Atoms representing the target system.

  * **opt_fn** (_Atoms - > Atoms_) – Optimization function to optimize the system after the deformation.

  * **show_logger** (_bool_) – Show log information.

  * **pbar** (_tqdm_ _or_ _None_) – Progress bar.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.



Returns
    

Resulting stress tensor.

Return type
    

[StressTensor](matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor "matlantis_features.features.elasticity.tensors.StressTensor")

[ __ previous matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms ](matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms.html "previous page") [ next matlantis_features.features.elasticity.tensors __](matlantis_features.features.elasticity.tensors.html "next page")
