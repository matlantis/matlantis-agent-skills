# matlantis_features.features.md.ase_integrators.ASEIntegrator#

_class _matlantis_features.features.md.ase_integrators.ASEIntegrator(_get_dynamics_fn : Callable[[Atoms], Dynamics]_, _timestep : float_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#ASEIntegrator)#
    

Bases: [`MDIntegratorBase`](matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase "matlantis_features.features.md.md_integrator_base.MDIntegratorBase")

Base class for the integrators using ASE backend.

Methods

`__init__`(get_dynamics_fn, timestep) | Initialize the integrator using ASE backend.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_get_dynamics_fn : Callable[[Atoms], Dynamics]_, _timestep : float_) → None[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#ASEIntegrator.__init__)#
    

Initialize the integrator using ASE backend.

Parameters
    

  * **get_dynamics_fn** (_CreateDynFunc_) – Function object for creating the Dynamics object of ASE.

  * **timestep** (_float_) – Integration time step in fs unit.




create_ase_dynamics(_atoms : Atoms_) → Dynamics[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#ASEIntegrator.create_ase_dynamics)#
    

Create the ASE’s Dynamics object.

Parameters
    

**atoms** (_Atoms_) – ASE’s Atoms object containing the system to simulate.

Returns
    

ASE’s Dynamics object.

Return type
    

Dynamics

set_temperature(_value : float_) → None[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#ASEIntegrator.set_temperature)#
    

Set the target temperature of integrator.

If this integrator does not support temperature control, this method raises the NotImplementedError exception.

Parameters
    

**value** (_float_) – Temperature in K unit.

[ __ previous matlantis_features.features.md.ase_integrators ](matlantis_features.features.md.ase_integrators.html "previous page") [ next matlantis_features.features.md.ase_integrators.AndersenIntegrator __](matlantis_features.features.md.ase_integrators.AndersenIntegrator.html "next page")
