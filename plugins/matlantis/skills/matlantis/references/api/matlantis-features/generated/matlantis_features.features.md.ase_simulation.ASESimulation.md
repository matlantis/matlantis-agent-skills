# matlantis_features.features.md.ase_simulation.ASESimulation#

_class _matlantis_features.features.md.ase_simulation.ASESimulation(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _estimator_fn : Optional[Callable[[], Estimator]] = None_)[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASESimulation)#
    

Bases: [`MDSimulationBase`](matlantis_features.features.md.md_simulation_base.MDSimulationBase.html#matlantis_features.features.md.md_simulation_base.MDSimulationBase "matlantis_features.features.md.md_simulation_base.MDSimulationBase")

MD Simulation object using ASE MD backend.

Methods

`__init__`(system, integrator[, estimator_fn]) | Create MD Simulation object using ASE MD backend.  
---|---  
`append_extension`(extension, interval) | Append extension object.  
`run`(n_run_steps) | Run the MD simulation.  
  
__init__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_, _estimator_fn : Optional[Callable[[], Estimator]] = None_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASESimulation.__init__)#
    

Create MD Simulation object using ASE MD backend.

Parameters
    

  * **system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system of the MD simulation. This object should be a subclass of ASEMDSystem

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation. This object should be a subclass of ASEIntegrator.

  * **estimator_fn** (_EstimatorFnType_ _or_ _None_ _,__optional_) – A factory method to create a custom estimator. Please refer [Customizing estimator used in matlantis-features](../tips/custom_estimator.html#custom-estimator) for detail.




append_extension(_extension : [MDExtensionBase](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")_, _interval : int_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASESimulation.append_extension)#
    

Append extension object.

Parameters
    

  * **extension** ([_MDExtensionBase_](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")) – Extension object appended to this simulation.

  * **interval** (_int_) – Timestep interval to invoke this extension object.




run(_n_run_steps : int_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASESimulation.run)#
    

Run the MD simulation.

Parameters
    

**n_run_steps** (_int_) – Number of simulation steps to be run.

[ __ previous matlantis_features.features.md.ase_simulation.ASECheckpoint ](matlantis_features.features.md.ase_simulation.ASECheckpoint.html "previous page") [ next matlantis_features.features.md.ase_simulation.ASETrajWriter __](matlantis_features.features.md.ase_simulation.ASETrajWriter.html "next page")
