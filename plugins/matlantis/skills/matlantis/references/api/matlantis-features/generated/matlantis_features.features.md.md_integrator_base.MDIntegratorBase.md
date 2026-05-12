# matlantis_features.features.md.md_integrator_base.MDIntegratorBase#

_class _matlantis_features.features.md.md_integrator_base.MDIntegratorBase(_timestep : float_)[[source]](../_modules/matlantis_features/features/md/md_integrator_base.html#MDIntegratorBase)#
    

Bases: `object`

Base class for md integrator.

Methods

`__init__`(timestep) | Initialize an instance.  
---|---  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_) → None[[source]](../_modules/matlantis_features/features/md/md_integrator_base.html#MDIntegratorBase.__init__)#
    

Initialize an instance.

Parameters
    

**timestep** (_float_) – Integration time step in fs unit.

set_temperature(_value : float_) → None[[source]](../_modules/matlantis_features/features/md/md_integrator_base.html#MDIntegratorBase.set_temperature)#
    

Set the target temperature of integrator.

If this integrator does not support temperature control, this method raises the NotImplementedError exception.

Parameters
    

**value** (_float_) – Temperature in K unit.

[ __ previous matlantis_features.features.md.md_integrator_base ](matlantis_features.features.md.md_integrator_base.html "previous page") [ next matlantis_features.features.md.md_simulation_base __](matlantis_features.features.md.md_simulation_base.html "next page")
