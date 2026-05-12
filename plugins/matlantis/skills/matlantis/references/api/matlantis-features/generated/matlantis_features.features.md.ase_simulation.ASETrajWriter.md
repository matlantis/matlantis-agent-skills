# matlantis_features.features.md.ase_simulation.ASETrajWriter#

_class _matlantis_features.features.md.ase_simulation.ASETrajWriter(_file_name : str_, _write_mode : str_, _props : List[str]_, _save_init_frame : bool = False_)[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASETrajWriter)#
    

Bases: [`MDExtensionBase`](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")

ASE’s trajectory writer extension.

Methods

`__init__`(file_name, write_mode, props[, ...]) | Create ASE's trajectory writer extension.  
---|---  
`__call__`(system, integrator) | Write the current trajectory frame.  
`close`() | Close the underlying trajectory stream.  
  
__init__(_file_name : str_, _write_mode : str_, _props : List[str]_, _save_init_frame : bool = False_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASETrajWriter.__init__)#
    

Create ASE’s trajectory writer extension.

Parameters
    

  * **file_name** (_str_) – Path name of the trajectory file.

  * **write_mode** (_str_) – Mode string (“w” or “a”, for write and append, respectively).

  * **props** (_list_ _[__str_ _]_) – List of property names to be written in the traj file.

  * **save_init_frame** (_bool_ _,__optional_) – If True, the initial structure of the simulation is also stored in the trajectory file. Defaults to False.




__call__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASETrajWriter.__call__)#
    

Write the current trajectory frame.

Parameters
    

  * **system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system of the MD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the MD simulation.




close() → None[[source]](../_modules/matlantis_features/features/md/ase_simulation.html#ASETrajWriter.close)#
    

Close the underlying trajectory stream.

[ __ previous matlantis_features.features.md.ase_simulation.ASESimulation ](matlantis_features.features.md.ase_simulation.ASESimulation.html "previous page") [ next matlantis_features.features.md.md __](matlantis_features.features.md.md.html "next page")
