# matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator#

_class _matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator(_timestep : float_, _temperature : float_, _taut : float = 500_, _fixcm : bool = True_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NVTBerendsenIntegrator)#
    

Bases: [`ASEIntegrator`](matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator "matlantis_features.features.md.ase_integrators.ASEIntegrator")

Berendsen thermostat integrator generating NVT ensemble, using the ASE backend.

Methods

`__init__`(timestep, temperature[, taut, fixcm]) | Create the Berendsen thermostat integrator.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_, _temperature : float_, _taut : float = 500_, _fixcm : bool = True_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NVTBerendsenIntegrator.__init__)#
    

Create the Berendsen thermostat integrator.

Parameters
    

  * **timestep** (_float_) – Integration time step in fs unit.

  * **temperature** (_float_) – Target temperature of the integrator in K unit.

  * **taut** (_float_ _,__optional_) – Time constant for Berendsen temperature coupling in fs. Defaults to 500 fs.

  * **fixcm** (_bool_ _,__optional_) – If True, the position and momentum of the center of mass is kept unperturbed. Defaults to True.




create_ase_dynamics(_atoms : Atoms_) → Dynamics#
    

Create the ASE’s Dynamics object.

Parameters
    

**atoms** (_Atoms_) – ASE’s Atoms object containing the system to simulate.

Returns
    

ASE’s Dynamics object.

Return type
    

Dynamics

set_temperature(_value : float_) → None#
    

Set the target temperature of integrator.

If this integrator does not support temperature control, this method raises the NotImplementedError exception.

Parameters
    

**value** (_float_) – Temperature in K unit.

[ __ previous matlantis_features.features.md.ase_integrators.NPTIntegrator ](matlantis_features.features.md.ase_integrators.NPTIntegrator.html "previous page") [ next matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator __](matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.html "next page")
