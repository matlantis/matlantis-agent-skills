# matlantis_features.features.md.md_extensions.DeformScheduler#

_class _matlantis_features.features.md.md_extensions.DeformScheduler(_final : ndarray_, _num_total_steps : int_)[[source]](../_modules/matlantis_features/features/md/md_extensions.html#DeformScheduler)#
    

Bases: [`MDExtensionBase`](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")

Linear cell shape scheduler extension.

Methods

`__init__`(final, num_total_steps) | Create linear cell shape scheduler.  
---|---  
`__call__`(system, integrator) | Change the cell of the system.  
  
__init__(_final : ndarray_, _num_total_steps : int_)[[source]](../_modules/matlantis_features/features/md/md_extensions.html#DeformScheduler.__init__)#
    

Create linear cell shape scheduler.

Parameters
    

  * **final** (_np.ndarray_) – The final cell shape after deformation.

  * **num_total_steps** (_int_) – Number of MD steps for the cell shape change procedure.




__call__(_system : [ASEMDSystem](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/md_extensions.html#DeformScheduler.__call__)#
    

Change the cell of the system.

Parameters
    

  * **system** ([_ASEMDSystem_](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")) – Target system of the MD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation.




[ __ previous matlantis_features.features.md.md_extensions ](matlantis_features.features.md.md_extensions.html "previous page") [ next matlantis_features.features.md.md_extensions.MDExtensionBase __](matlantis_features.features.md.md_extensions.MDExtensionBase.html "next page")
