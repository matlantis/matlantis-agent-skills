# matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension#

_class _matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension(_rnemd_type : str_, _n_slab : int_)[[source]](../_modules/matlantis_features/features/md/post_features/rnemd_extension.html#RNEMDExtension)#
    

Bases: [`MDExtensionBase`](matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase "matlantis_features.features.md.md_extensions.MDExtensionBase")

Extension class for the reverse non-equilibrium MD simulation.

Methods

`__init__`(rnemd_type, n_slab) | Initialize an instance.  
---|---  
`__call__`(system, integrator) | Perform the exchange process of the rNEMD simulation.  
  
__init__(_rnemd_type : str_, _n_slab : int_)[[source]](../_modules/matlantis_features/features/md/post_features/rnemd_extension.html#RNEMDExtension.__init__)#
    

Initialize an instance.

Parameters
    

  * **rnemd_type** (_str_) – Name of the rNEMD calculation type

  * **"thermal_conductivity"****)****.** (_(__"viscosity" or_) – 

  * **n_slab** (_int_) – Number of the slabs.




__call__(_system : [MDSystemBase](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")_, _integrator : [MDIntegratorBase](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")_) → None[[source]](../_modules/matlantis_features/features/md/post_features/rnemd_extension.html#RNEMDExtension.__call__)#
    

Perform the exchange process of the rNEMD simulation.

Parameters
    

  * **system** ([_MDSystemBase_](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")) – Target system of the rNEMD simulation.

  * **integrator** ([_MDIntegratorBase_](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")) – Integrator used for the rNEMD simulation.




[ __ previous matlantis_features.features.md.post_features.rnemd_extension ](matlantis_features.features.md.post_features.rnemd_extension.html "previous page") [ next matlantis_features.features.md.post_features.rnemd_extension.fit_slope __](matlantis_features.features.md.post_features.rnemd_extension.fit_slope.html "next page")
