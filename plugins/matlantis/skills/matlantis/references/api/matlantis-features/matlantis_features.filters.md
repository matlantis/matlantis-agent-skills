# matlantis_features.filters#

_class _matlantis_features.filters.ExpCellASEFilter(_mask : Optional[List[bool]] = None_, _cell_factor : Optional[float] = None_, _hydrostatic_strain : bool = False_, _constant_volume : bool = False_, _scalar_pressure : float = 0.0_)[[source]](_modules/matlantis_features/filters.html#ExpCellASEFilter)#
    

Bases: `MatlantisASEFilter`

The filter that uses the ASE ExpCellASEFilter implemmentation.

__init__(_mask : Optional[List[bool]] = None_, _cell_factor : Optional[float] = None_, _hydrostatic_strain : bool = False_, _constant_volume : bool = False_, _scalar_pressure : float = 0.0_)[[source]](_modules/matlantis_features/filters.html#ExpCellASEFilter.__init__)#
    

Initialize an instance.

Parameters
    

  * **mask** (_list_ _[__bool_ _] or_ _None_ _,__optional_) – A list of booleans indicating which of the six independent components of the strain are relaxed. If None, all components are relaxed. Defaults to None.

  * **cell_factor** (_float_ _or_ _None_ _,__optional_) – Factor by which deformation gradient is multiplied to put it on the same scale as the positions when assembling the combined position/cell vector. The stress contribution to the forces is scaled down by the same factor. If None, a default value which is equal to the number of atoms will be used. Defaults to None.

  * **hydrostatic_strain** (_bool_ _,__optional_) – Constrain the cell by only allowing hydrostatic deformation. Defaults to False.

  * **constant_volume** (_bool_ _,__optional_) – Relax cell while keep the cell volume. Note: this only approximately conserves the volume and breaks energy/force consistency so can only be used with optimizers that do require do a line minimisation (e.g. FIRE). Defaults to False.

  * **scalar_pressure** (_float_ _,__optional_) – Applied pressure (eV/Angstrom^3). As above, this breaks energy/force consistency. Defaults to 0.0.




to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/filters.html#ExpCellASEFilter.to_dict)#
    

Dictionary representation of the ExpCellASEFilter.

Returns
    

A dict containing a serialized ExpCellASEFilter.

Return type
    

dict[str, Any]

_class _matlantis_features.filters.IsotropicASEFilter[[source]](_modules/matlantis_features/filters.html#IsotropicASEFilter)#
    

Bases: `MatlantisASEFilter`

The filter that uses the IsotropicFilter implemmentation.

__init__() → None[[source]](_modules/matlantis_features/filters.html#IsotropicASEFilter.__init__)#
    

Initialize an instance.

to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/filters.html#IsotropicASEFilter.to_dict)#
    

Dictionary representation of the IsotropicASEFilter.

Returns
    

A dict containing a serialized IsotropicASEFilter.

Return type
    

dict[str, Any]

_class _matlantis_features.filters.IsotropicFilter(_atoms : Atoms_)[[source]](_modules/matlantis_features/filters.html#IsotropicFilter)#
    

Bases: `Filter`

The filter to keep the ratio of cell lengths during MD or optimization.

__init__(_atoms : Atoms_) → None[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.__init__)#
    

Initialize an instance.

Parameters
    

**atoms** (_ASEAtoms_) – Filter atoms.

get_global_number_of_atoms() → int[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.get_global_number_of_atoms)#
    

Returns the global number of atoms in a distributed-atoms parallel simulation.

Returns
    

The global number of atoms.

Return type
    

int

get_kinetic_energy() → float[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.get_kinetic_energy)#
    

Get the kinetic energy.

Returns
    

The kinetic energy.

Return type
    

float

get_stress(_** kwargs: Any_) → ndarray[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.get_stress)#
    

Calculate stress tensor.

Returns
    

Stress tensor after applying filter.

Return type
    

np.ndarray

get_temperature() → float[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.get_temperature)#
    

Get the temperature in Kelvin.

Returns
    

The temperature.

Return type
    

float

get_total_energy() → float[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.get_total_energy)#
    

Get the total energy.

Returns
    

The total energy.

Return type
    

float

set_cell(_cell : ndarray_, _scale_atoms : bool = False_, _apply_constraint : bool = True_) → None[[source]](_modules/matlantis_features/filters.html#IsotropicFilter.set_cell)#
    

Set unit cell vectors.

Parameters
    

  * **cell** (_np.ndarray_) – Unit cell.

  * **scale_atoms** (_bool_ _,__optional_) – Whether to fix atomic positions or move atoms with the unit cell. Defaults to False.

  * **apply_constraint** (_bool_ _,__optional_) – Whether to apply constraints to the given cell. Defaults to True.




_class _matlantis_features.filters.MatlantisASEFilter(_get_filter : Callable[[Union[Atoms, [MatlantisAtoms](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]], Filter]_)[[source]](_modules/matlantis_features/filters.html#MatlantisASEFilter)#
    

Bases: `MatlantisFilter`

MatlantisFilter class will be used in the optimization feature.

__call__(_atoms : Union[Atoms, [MatlantisAtoms](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_) → Filter[[source]](_modules/matlantis_features/filters.html#MatlantisASEFilter.__call__)#
    

Call function.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – The system to be optimized. The ASEAtoms and MatlantisAtoms are also supported.

Returns
    

The filter to be used in the optimization.

Return type
    

Filter

__init__(_get_filter : Callable[[Union[Atoms, [MatlantisAtoms](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]], Filter]_)[[source]](_modules/matlantis_features/filters.html#MatlantisASEFilter.__init__)#
    

Initialize an instance.

Parameters
    

**get_filter** (_ASEAtoms_ _or_ _MatlantisAtoms - > Filter_) – A function to get a filter.

_class _matlantis_features.filters.MatlantisFilter[[source]](_modules/matlantis_features/filters.html#MatlantisFilter)#
    

Bases: `object`

Template base class for Matlantis Filter.

__call__(_atoms : Union[Atoms, [MatlantisAtoms](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")]_) → Any[[source]](_modules/matlantis_features/filters.html#MatlantisFilter.__call__)#
    

Call function.

Parameters
    

**atoms** (_ASEAtoms_ _or_[ _MatlantisAtoms_](matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms")) – The system to add filter to.

Returns
    

The system with filter which to be used in the optimization.

Return type
    

Any

_classmethod _from_dict(_res : Dict[str, Any]_) → MatlantisFilter[[source]](_modules/matlantis_features/filters.html#MatlantisFilter.from_dict)#
    

Construct a MatlantisFilter object from serialized dict.

Parameters
    

**res** (_dict_ _[__str_ _,__Any_ _]_) – A dict containing a serialized MatlantisFilter from to_dict().

Returns
    

The deserialized object from provided dict.

Return type
    

MatlantisFilter

to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/filters.html#MatlantisFilter.to_dict)#
    

Dictionary representation of the MatlantisFilter.

Returns
    

A dict containing a serialized MatlantisFilter.

Return type
    

dict[str, Any]

_class _matlantis_features.filters.StrainASEFilter(_mask : Optional[List[bool]] = None_, _include_ideal_gas : bool = False_)[[source]](_modules/matlantis_features/filters.html#StrainASEFilter)#
    

Bases: `MatlantisASEFilter`

The filter that uses the ASE StrainASEFilter implemmentation.

__init__(_mask : Optional[List[bool]] = None_, _include_ideal_gas : bool = False_)[[source]](_modules/matlantis_features/filters.html#StrainASEFilter.__init__)#
    

Initialize an instance.

Parameters
    

  * **mask** (_list_ _[__bool_ _] or_ _None_ _,__optional_) – A list of booleans indicating which of the six independent components of the strain are relaxed. If None, all components are relaxed. Defaults to None.

  * **include_ideal_gas** (_bool_ _,__optional_) – The ideal gas contribution to the stresses is added. Defaults to False.




to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/filters.html#StrainASEFilter.to_dict)#
    

Dictionary representation of the StrainASEFilter.

Returns
    

A dict containing a serialized StrainASEFilter.

Return type
    

dict[str, Any]

_class _matlantis_features.filters.UnitCellASEFilter(_mask : Optional[List[bool]] = None_, _cell_factor : Optional[float] = None_, _hydrostatic_strain : bool = False_, _constant_volume : bool = False_, _scalar_pressure : float = 0.0_)[[source]](_modules/matlantis_features/filters.html#UnitCellASEFilter)#
    

Bases: `MatlantisASEFilter`

The filter that uses the ASE UnitCellFilter implemmentation.

__init__(_mask : Optional[List[bool]] = None_, _cell_factor : Optional[float] = None_, _hydrostatic_strain : bool = False_, _constant_volume : bool = False_, _scalar_pressure : float = 0.0_)[[source]](_modules/matlantis_features/filters.html#UnitCellASEFilter.__init__)#
    

Initialize an instance.

Parameters
    

  * **mask** (_list_ _[__bool_ _] or_ _None_ _,__optional_) – A list of booleans indicating which of the six independent components of the strain are relaxed. Defaults to None.

  * **cell_factor** (_float_ _or_ _None_ _,__optional_) – Factor by which deformation gradient is multiplied to put it on the same scale as the positions when assembling the combined position/cell vector. The stress contribution to the forces is scaled down by the same factor. If None, a default value which is equal to the number of atoms will be used. Defaults to None.

  * **hydrostatic_strain** (_bool_ _,__optional_) – Constrain the cell by only allowing hydrostatic deformation. Defaults to False.

  * **constant_volume** (_bool_ _,__optional_) – Relax cell while keep the cell volume. Note: this only approximately conserves the volume and breaks energy/force consistency so can only be used with optimizers that do require do a line minimisation (e.g. FIRE). Defaults to False.

  * **scalar_pressure** (_float_ _,__optional_) – Applied pressure (eV/Angstrom^3). As above, this breaks energy/force consistency. Defaults to 0.0.




to_dict() → Dict[str, Any][[source]](_modules/matlantis_features/filters.html#UnitCellASEFilter.to_dict)#
    

Dictionary representation of the UnitCellASEFilter.

Returns
    

A dict containing a serialized UnitCellASEFilter.

Return type
    

dict[str, Any]

[ __ previous matlantis_features.atoms ](matlantis_features.atoms.html "previous page") [ next matlantis_features.ase_ext __](matlantis_features.ase_ext.html "next page")
