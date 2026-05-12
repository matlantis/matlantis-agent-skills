# matlantis_features.features.md.ase_integrators.NPTIntegrator#

_class _matlantis_features.features.md.ase_integrators.NPTIntegrator(_timestep : float_, _temperature : float_, _pressure : float_, _ttime : float_, _pfactor : float_, _mask : Optional[ndarray] = None_, _hydrostatic_strain : bool = False_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NPTIntegrator)#
    

Bases: [`ASEIntegrator`](matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator "matlantis_features.features.md.ase_integrators.ASEIntegrator")

Nose-Hoover/Parrinello-Rahman dynamics integrator generating NPT ensemble, using the ASE backend.

Methods

`__init__`(timestep, temperature, pressure, ...) | Create the Nose-Hoover/Parrinello-Rahman dynamics integrator.  
---|---  
`create_ase_dynamics`(atoms) | Create the ASE's Dynamics object.  
`set_temperature`(value) | Set the target temperature of integrator.  
  
__init__(_timestep : float_, _temperature : float_, _pressure : float_, _ttime : float_, _pfactor : float_, _mask : Optional[ndarray] = None_, _hydrostatic_strain : bool = False_)[[source]](../_modules/matlantis_features/features/md/ase_integrators.html#NPTIntegrator.__init__)#
    

Create the Nose-Hoover/Parrinello-Rahman dynamics integrator.

Parameters
    

  * **timestep** (_float_) – Integration time step in fs unit.

  * **temperature** (_float_) – Target temperature of the integrator in K unit.

  * **pressure** (_float_) – Target pressure in eV/Angstrom^3 unit.

  * **ttime** (_float_) – Characteristic timescale of the thermostat, in fs units.

  * **pfactor** (_float_) – A constant in the barostat differential equation. If a characteristic barostat timescale of ptime is desired, set pfactor to ptime^2 * B (where ptime is in units matching eV, Angstrom, u; and B is the Bulk Modulus, given in eV/Angstrom^3).

  * **mask** (_np.ndarray_ _or_ _None_ _,__optional_) – A 3x3 array indicating if the system can change size along the three Cartesian axes. Defaults to None.

  * **hydrostatic_strain** (_bool_ _,__optional_) – Constrain the cell by only allowing hydrostatic deformation. The virial tensor is replaced by np.diag([np.trace(virial)]*3). Defaults to False.




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

[ __ previous matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator ](matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator.html "previous page") [ next matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator __](matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator.html "next page")
