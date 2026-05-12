# matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator#

_class _matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator(_timestep : float_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#VelocityVerletIntegrator)#
    

Bases: [`ASEIntegrator`](matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator "matlantis_features.features.md.ase_integrators.ASEIntegrator")

Velocity Verlet integrator using ASE backend.

Methods

`__init__`(timestep) | Create a Velocity Verlet integrator.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_) → None[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#VelocityVerletIntegrator.__init__)#
    

Create a Velocity Verlet integrator.

Parameters
    

**timestep** (_float_) – Integration time step in fs unit.

create_ase_dynamics(_atoms : Atoms_) → Dynamics#
    

Create the ASE’s Dynamics object.

Parameters
    

**atoms** (_Atoms_) – ASE’s Atoms object containing the system to simulate.

Returns
    

ASE’s Dynamics object.

Return type
    

Dynamics

set_temperature(_value : float_) → None[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#VelocityVerletIntegrator.set_temperature)#
    

Set the target temperature of integrator.

This integrator does not support temperature control, so this method always raises the NotImplementedError exception.

Parameters
    

**value** (_float_) – Temperature in K unit.

[ __ previous matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator ](matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator.html "previous page") [ next matlantis_features.features.md.ase_md_system __](matlantis_features.features.md.ase_md_system.html "next page")
