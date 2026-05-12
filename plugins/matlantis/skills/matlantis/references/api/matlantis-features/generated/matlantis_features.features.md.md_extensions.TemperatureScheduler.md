# matlantis_features.features.md.md_extensions.TemperatureScheduler#

_class _matlantis_features.features.md.md_extensions.TemperatureScheduler(_start_value : float_, _end_value : float_, _num_total_steps : int_)[[source]](../_modules/matlantis_features/features/md/md_extensions.html#TemperatureScheduler)#
    

Bases: [`MDExtensionBase`](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")

Linear temperature scheduler exteisnon.

This extension linearly change the temperature of the system from the start to the end temperature, during the given time steps. After the given time step, the temperature is kept to the end temperature.

Methods

`__init__`(start_value, end_value, num_total_steps) | Create linear temperature scheduler exteisnon.  
---|---  
`__call__`(system, integrator) | Change the temperature of the system.  
  
__init__(_start_value : float_, _end_value : float_, _num_total_steps : int_)[[source]](../_modules/matlantis_features/features/md/md_extensions.html#TemperatureScheduler.__init__)#
    

Create linear temperature scheduler exteisnon.

Parameters
    

  * **start_value** (_float_) – Start temperature in K unit.

  * **end_value** (_float_) – End temperature in K unit.

  * **num_total_steps** (_int_) – Number of steps for the linear temperature change.




__call__(_system : [ASEMDSystem](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/md_extensions.html#TemperatureScheduler.__call__)#
    

Change the temperature of the system.

Parameters
    

  * **system** ([_ASEMDSystem_](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")) – Target system of the MD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation.




[ __ previous matlantis_features.features.md.md_extensions.MDExtensionBase ](matlantis_features.features.md.md_extensions.MDExtensionBase.html "previous page") [ next matlantis_features.features.md.md_integrator_base __](matlantis_features.features.md.md_integrator_base.html "next page")
