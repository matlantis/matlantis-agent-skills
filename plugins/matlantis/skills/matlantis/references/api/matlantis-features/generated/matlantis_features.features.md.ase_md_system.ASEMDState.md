# matlantis_features.features.md.ase_md_system.ASEMDState#

_class _matlantis_features.features.md.ase_md_system.ASEMDState(_system_state : Dict[str, Any]_, _integrator_state : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDState)#
    

Bases: `object`

Serializable MD simulation’s state class.

Methods

`__init__`(system_state, integrator_state) | Create the MD state object.  
---|---  
`from_file`(file_name) | Create the MD state object from file.  
`restore_to`(system, integrator) | Restore the state data stored in this object to the specified system and integrator.  
  
__init__(_system_state : Dict[str, Any]_, _integrator_state : Dict[str, Any]_)[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDState.__init__)#
    

Create the MD state object.

Parameters
    

  * **system_state** (_MDStateDict_) – System’s state object.

  * **integrator_state** (_MDStateDict_) – Integrator’s state object.




_static _from_file(_file_name : str_) → ASEMDState[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDState.from_file)#
    

Create the MD state object from file.

Parameters
    

**file_name** (_str_) – Path name of the checkpoint file to load.

Returns
    

Created the MD state object.

Return type
    

ASEMDState

restore_to(_system : [ASEMDSystem](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDState.restore_to)#
    

Restore the state data stored in this object to the specified system and integrator.

Parameters
    

  * **system** ([_ASEMDSystem_](matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem "matlantis_features.features.md.ase_md_system.ASEMDSystem")) – Target system object.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Target integrator object.




[ __ previous matlantis_features.features.md.ase_md_system ](matlantis_features.features.md.ase_md_system.html "previous page") [ next matlantis_features.features.md.ase_md_system.ASEMDSystem __](matlantis_features.features.md.ase_md_system.ASEMDSystem.html "next page")
