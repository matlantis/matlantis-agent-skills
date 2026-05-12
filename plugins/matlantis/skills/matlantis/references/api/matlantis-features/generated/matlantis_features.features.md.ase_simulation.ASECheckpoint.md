# matlantis_features.features.md.ase_simulation.ASECheckpoint#

_class _matlantis_features.features.md.ase_simulation.ASECheckpoint(_file_name : str_)[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASECheckpoint)#
    

Bases: [`MDExtensionBase`](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")

Checkpoint extension for resuming the ASE MD simulation.

Methods

`__init__`(file_name) | Create checkpoint extension for resuming the ASE MD simulation.  
---|---  
`__call__`(system, integrator) | Create the checkpoint file of the current simulation state.  
  
__init__(_file_name : str_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASECheckpoint.__init__)#
    

Create checkpoint extension for resuming the ASE MD simulation.

Parameters
    

**file_name** (_str_) – Path of the checkpoint file to be created.

__call__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASECheckpoint.__call__)#
    

Create the checkpoint file of the current simulation state.

Parameters
    

  * **system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system of the MD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation.




[ __ previous matlantis_features.features.md.ase_simulation ](matlantis_features.features.md.ase_simulation.html "previous page") [ next matlantis_features.features.md.ase_simulation.ASESimulation __](matlantis_features.features.md.ase_simulation.ASESimulation.html "next page")
