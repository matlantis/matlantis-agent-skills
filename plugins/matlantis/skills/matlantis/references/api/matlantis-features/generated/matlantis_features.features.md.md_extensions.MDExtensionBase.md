# matlantis_features.features.md.md_extensions.MDExtensionBase#

_class _matlantis_features.features.md.md_extensions.MDExtensionBase[[source]](../_modules/matlantis_features/features/md/md_extensions.html#MDExtensionBase)#
    

Bases: `object`

Base class for the MD extension.

Methods

`__init__`() |   
---|---  
`__call__`(system, integrator) | Call the MD extension.  
  
__init__()#
    

__call__(_system : [ASEMDSystem](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/md_extensions.html#MDExtensionBase.__call__)#
    

Call the MD extension.

Parameters
    

  * **system** ([_ASEMDSystem_](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")) – Target system of the MD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation.




[ __ previous matlantis_features.features.md.md_extensions.DeformScheduler ](matlantis_features.features.md.md_extensions.DeformScheduler.html "previous page") [ next matlantis_features.features.md.md_extensions.TemperatureScheduler __](matlantis_features.features.md.md_extensions.TemperatureScheduler.html "next page")
