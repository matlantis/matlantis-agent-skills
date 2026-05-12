# matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator#

_class _matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator(_timestep : float_, _temperature : float_, _pressure : float_, _taut : float = 500_, _taup : float = 1000_, _compressibility : float = 67.2_, _fixcm : bool = True_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NPTBerendsenIntegrator)#
    

Bases: [`ASEIntegrator`](matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator "matlantis_features.features.md.ase_integrators.ASEIntegrator")

Berendsen thermo/barostat integrator generating an NPT ensemble via ASE.

Relax both temperature and isotropic pressure.

Methods

`__init__`(timestep, temperature, pressure[, ...]) | Create Berendsen thermo/barostat integrator.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_, _temperature : float_, _pressure : float_, _taut : float = 500_, _taup : float = 1000_, _compressibility : float = 67.2_, _fixcm : bool = True_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NPTBerendsenIntegrator.__init__)#
    

Create Berendsen thermo/barostat integrator.

Parameters
    

  * **timestep** (_float_) – Integration time step in fs.

  * **temperature** (_float_) – Target temperature in K.

  * **pressure** (_float_) – Target isotropic pressure in eV/Angstrom^3.

  * **taut** (_float_ _,__optional_) – Temperature coupling time constant (fs). Smaller values enforce temperature more strongly; defaults to 500 fs.

  * **taup** (_float_ _,__optional_) – Time constant for Berendsen pressure coupling. Defaults to 1 ps.

  * **compressibility** (_float_ _,__optional_) – The compressibility of the material (Angstrom^3/eV). Defaults to 67.2.

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

[ __ previous matlantis_features.features.md.ase_integrators.LangevinIntegrator ](matlantis_features.features.md.ase_integrators.LangevinIntegrator.html "previous page") [ next matlantis_features.features.md.ase_integrators.NPTIntegrator __](matlantis_features.features.md.ase_integrators.NPTIntegrator.html "next page")
