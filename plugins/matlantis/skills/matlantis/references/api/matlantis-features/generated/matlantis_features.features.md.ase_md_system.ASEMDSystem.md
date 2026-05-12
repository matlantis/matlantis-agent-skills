# matlantis_features.features.md.ase_md_system.ASEMDSystem#

_class _matlantis_features.features.md.ase_md_system.ASEMDSystem(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _step : int = 0_, _time : float = 0.0_)[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDSystem)#
    

Bases: [`MDSystemBase`](matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase "matlantis_features.features.md.md_system_base.MDSystemBase")

MD System using ASE backend.

Methods

`__init__`(atoms[, step, time]) | Create MD System using ASE backend.  
---|---  
`init_temperature`(temperature[, stationary, ...]) | Initialize the velocities by the Maxwell-Boltzmann distribution of the specified temperature.  
  
__init__(_atoms : Union[Atoms, [MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_, _step : int = 0_, _time : float = 0.0_)[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDSystem.__init__)#
    

Create MD System using ASE backend.

Parameters
    

  * **atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – ASE’s Atoms object representing the simulation system.

  * **step** (_int_ _,__optional_) – How many MD steps have been run before the task. Defaults to 0.

  * **time** (_float_ _,__optional_) – How long time (in fs) has it been simulated before the task. Defaults to 0.0.




init_temperature(_temperature : float_, _stationary : bool = False_, _zero_rotation : bool = False_, _rng : Optional[RandomState] = None_) → None[[source]](../_modules/matlantis_features/features/md/ase_md_system.html#ASEMDSystem.init_temperature)#
    

Initialize the velocities by the Maxwell-Boltzmann distribution of the specified temperature.

Parameters
    

  * **temperature** (_float_) – Target temperature in K unit.

  * **stationary** (_bool_ _,__optional_) – If True, set the center-of-mass momentum to zero. Defaults to False.

  * **zero_rotation** (_bool_ _,__optional_) – If True, set the total angular momentum to zero by counteracting rigid rotations. Defaults to False.

  * **rng** (_RandomState_ _or_ _None_ _,__optional_) – Random number generator to generate the distribution. If None, the default generator of ASE implementation is used. Defaults to None.




[ __ previous matlantis_features.features.md.ase_md_system.ASEMDState ](matlantis_features.features.md.ase_md_system.ASEMDState.html "previous page") [ next matlantis_features.features.md.ase_simulation __](matlantis_features.features.md.ase_simulation.html "next page")
