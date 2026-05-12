# matlantis_features.features.md.ase_integrators.AndersenIntegrator#

_class _matlantis_features.features.md.ase_integrators.AndersenIntegrator(_timestep : float_, _temperature : float_, _andersen_prob : float = 0.01_, _fixcm : bool = True_, _random_seed : int = 0_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#AndersenIntegrator)#
    

Bases: [`ASEIntegrator`](matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator "matlantis_features.features.md.ase_integrators.ASEIntegrator")

Andersen thermostat integrator generating NVT ensemble, using the ASE backend.

Methods

`__init__`(timestep, temperature[, ...]) | Create the Andersen integrator.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_, _temperature : float_, _andersen_prob : float = 0.01_, _fixcm : bool = True_, _random_seed : int = 0_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#AndersenIntegrator.__init__)#
    

Create the Andersen integrator.

Parameters
    

  * **timestep** (_float_) – Integration time step in fs unit.

  * **temperature** (_float_) – Target temperature of the integrator in K unit.

  * **andersen_prob** (_float_ _,__optional_) – A random collision probability (0-1). with this probability, atoms get assigned random velocity components. Defaults to 0.01.

  * **fixcm** (_bool_ _,__optional_) – If True, the position and momentum of the center of mass is kept unperturbed. Defaults to True.

  * **random_seed** (_int_ _,__optional_) – Random number generator’s seed used for the simulation. Defaults to 0.




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

[ __ previous matlantis_features.features.md.ase_integrators.ASEIntegrator ](matlantis_features.features.md.ase_integrators.ASEIntegrator.html "previous page") [ next matlantis_features.features.md.ase_integrators.LangevinIntegrator __](matlantis_features.features.md.ase_integrators.LangevinIntegrator.html "next page")
