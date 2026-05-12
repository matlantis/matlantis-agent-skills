# Source code for matlantis_features.filters
    
    
    from typing import Any, Callable, Dict, List, Optional, Union
    
    import ase
    import numpy as np
    from packaging import version
    
    if version.parse(ase.__version__) >= version.parse("3.23.0"):
        from ase.filters import ExpCellFilter, Filter, StrainFilter, UnitCellFilter
    else:
        from ase.constraints import ExpCellFilter, Filter, StrainFilter, UnitCellFilter
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.utils.save_util import from_dict
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter)class MatlantisFilter(object):
        """Template base class for Matlantis Filter."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter.__call__)    def __call__(self, atoms: Union[ASEAtoms, MatlantisAtoms]) -> Any:
            """Call function.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms): The system to add filter to.
            Returns:
                Any : The system with filter which to be used in the optimization.
            """
            raise NotImplementedError
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the MatlantisFilter.
    
            Returns:
                dict[str, Any] : A dict containing a serialized MatlantisFilter.
            """
            raise NotImplementedError
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisFilter.from_dict)    @classmethod
        def from_dict(cls, res: Dict[str, Any]) -> "MatlantisFilter":
            """Construct a MatlantisFilter object from serialized dict.
    
            Args:
                res (dict[str, Any]): A dict containing a serialized MatlantisFilter from to_dict().
            Returns:
                MatlantisFilter : The deserialized object from provided dict.
            """
            obj = from_dict(res)
            if not isinstance(obj, cls):
                raise TypeError(
                    f"The class {res['cls_path']} specified in 'cls_path' does not match the instantiated class ({type(obj).__name__})."
                )
            return obj
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisASEFilter)class MatlantisASEFilter(MatlantisFilter):
        """MatlantisFilter class will be used in the optimization feature."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisASEFilter.__init__)    def __init__(self, get_filter: Callable[[Union[ASEAtoms, MatlantisAtoms]], Filter]):
            """Initialize an instance.
    
            Args:
                get_filter (ASEAtoms or MatlantisAtoms -> Filter):  A function to get a filter.
            """
            self._get_filter = get_filter
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.MatlantisASEFilter.__call__)    def __call__(self, atoms: Union[ASEAtoms, MatlantisAtoms]) -> Filter:
            """Call function.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms):
                    The system to be optimized. The ASEAtoms and MatlantisAtoms are also supported.
            Returns:
                Filter : The filter to be used in the optimization.
            """
            ase_filter = self._get_filter(atoms)
            return ase_filter
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.UnitCellASEFilter)class UnitCellASEFilter(MatlantisASEFilter):
        """The filter that uses the ASE UnitCellFilter implemmentation."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.UnitCellASEFilter.__init__)    def __init__(
            self,
            mask: Optional[List[bool]] = None,
            cell_factor: Optional[float] = None,
            hydrostatic_strain: bool = False,
            constant_volume: bool = False,
            scalar_pressure: float = 0.0,
        ):
            """Initialize an instance.
    
            Args:
                mask (list[bool] or None, optional):
                    A list of booleans indicating which of the six independent components of the
                    strain are relaxed. Defaults to None.
                cell_factor (float or None, optional):
                    Factor by which deformation gradient is multiplied to put it on the same scale as
                    the positions when assembling the combined position/cell vector. The stress
                    contribution to the forces is scaled down by the same factor. If None, a default
                    value which is equal to the number of atoms will be used. Defaults to None.
                hydrostatic_strain (bool, optional):
                    Constrain the cell by only allowing hydrostatic deformation. Defaults to False.
                constant_volume (bool, optional):
                    Relax cell while keep the cell volume.
                    Note: this only approximately conserves the volume and breaks energy/force
                    consistency so can only be used with optimizers that do require do a line
                    minimisation (e.g. FIRE). Defaults to False.
                scalar_pressure (float, optional):
                    Applied pressure (eV/Angstrom^3). As above, this breaks
                    energy/force consistency. Defaults to 0.0.
            """
            self.mask = mask
            self.cell_factor = cell_factor
            self.hydrostatic_strain = hydrostatic_strain
            self.constant_volume = constant_volume
            self.scalar_pressure = scalar_pressure
    
            def get_filter(atoms: Union[ASEAtoms, MatlantisAtoms]) -> Filter:
                if isinstance(atoms, MatlantisAtoms):
                    system = atoms.ase_atoms
                else:
                    system = atoms
    
                return UnitCellFilter(
                    system,
                    mask=mask,
                    cell_factor=cell_factor,
                    hydrostatic_strain=hydrostatic_strain,
                    constant_volume=constant_volume,
                    scalar_pressure=scalar_pressure,
                )
    
            super(UnitCellASEFilter, self).__init__(get_filter=get_filter)
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.UnitCellASEFilter.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the UnitCellASEFilter.
    
            Returns:
                dict[str, Any] : A dict containing a serialized UnitCellASEFilter.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "mask": self.mask,
                "cell_factor": self.cell_factor,
                "hydrostatic_strain": self.hydrostatic_strain,
                "constant_volume": self.constant_volume,
                "scalar_pressure": self.scalar_pressure,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.StrainASEFilter)class StrainASEFilter(MatlantisASEFilter):
        """The filter that uses the ASE StrainASEFilter implemmentation."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.StrainASEFilter.__init__)    def __init__(
            self,
            mask: Optional[List[bool]] = None,
            include_ideal_gas: bool = False,
        ):
            """Initialize an instance.
    
            Args:
                mask (list[bool] or None, optional):
                    A list of booleans indicating which of the six independent components of the
                    strain are relaxed. If None, all components are relaxed. Defaults to None.
                include_ideal_gas (bool, optional):
                    The ideal gas contribution to the stresses is added. Defaults to False.
            """
            self.mask = mask
            self.include_ideal_gas = include_ideal_gas
    
            def get_filter(atoms: Union[ASEAtoms, MatlantisAtoms]) -> Filter:
                if isinstance(atoms, MatlantisAtoms):
                    system = atoms.ase_atoms
                else:
                    system = atoms
    
                return StrainFilter(
                    system,
                    mask=mask,
                    include_ideal_gas=include_ideal_gas,
                )
    
            super(StrainASEFilter, self).__init__(get_filter=get_filter)
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.StrainASEFilter.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the StrainASEFilter.
    
            Returns:
                dict[str, Any] : A dict containing a serialized StrainASEFilter.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "mask": self.mask,
                "include_ideal_gas": self.include_ideal_gas,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.ExpCellASEFilter)class ExpCellASEFilter(MatlantisASEFilter):
        """The filter that uses the ASE ExpCellASEFilter implemmentation."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.ExpCellASEFilter.__init__)    def __init__(
            self,
            mask: Optional[List[bool]] = None,
            cell_factor: Optional[float] = None,
            hydrostatic_strain: bool = False,
            constant_volume: bool = False,
            scalar_pressure: float = 0.0,
        ):
            """Initialize an instance.
    
            Args:
                mask (list[bool] or None, optional):
                    A list of booleans indicating which of the six independent components of the
                    strain are relaxed. If None, all components are relaxed. Defaults to None.
                cell_factor (float or None, optional):
                    Factor by which deformation gradient is multiplied to put it on the same scale as
                    the positions when assembling the combined position/cell vector. The stress
                    contribution to the forces is scaled down by the same factor. If None, a default
                    value which is equal to the number of atoms will be used. Defaults to None.
                hydrostatic_strain (bool, optional):
                    Constrain the cell by only allowing hydrostatic deformation. Defaults to False.
                constant_volume (bool, optional):
                    Relax cell while keep the cell volume.
                    Note: this only approximately conserves the volume and breaks energy/force
                    consistency so can only be used with optimizers that do require do a line
                    minimisation (e.g. FIRE). Defaults to False.
                scalar_pressure (float, optional):
                    Applied pressure (eV/Angstrom^3). As above, this breaks
                    energy/force consistency. Defaults to 0.0.
            """
            self.mask = mask
            self.cell_factor = cell_factor
            self.hydrostatic_strain = hydrostatic_strain
            self.constant_volume = constant_volume
            self.scalar_pressure = scalar_pressure
    
            def get_filter(atoms: Union[ASEAtoms, MatlantisAtoms]) -> Filter:
                if isinstance(atoms, MatlantisAtoms):
                    system = atoms.ase_atoms
                else:
                    system = atoms
    
                return ExpCellFilter(
                    system,
                    mask=mask,
                    cell_factor=cell_factor,
                    hydrostatic_strain=hydrostatic_strain,
                    constant_volume=constant_volume,
                    scalar_pressure=scalar_pressure,
                )
    
            super(ExpCellASEFilter, self).__init__(get_filter=get_filter)
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.ExpCellASEFilter.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the ExpCellASEFilter.
    
            Returns:
                dict[str, Any] : A dict containing a serialized ExpCellASEFilter.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "mask": self.mask,
                "cell_factor": self.cell_factor,
                "hydrostatic_strain": self.hydrostatic_strain,
                "constant_volume": self.constant_volume,
                "scalar_pressure": self.scalar_pressure,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter)class IsotropicFilter(Filter):  # type: ignore
        """The filter to keep the ratio of cell lengths during MD or optimization."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.__init__)    def __init__(self, atoms: ASEAtoms) -> None:
            """Initialize an instance.
    
            Args:
                atoms (ASEAtoms): Filter atoms.
            """
            super(IsotropicFilter, self).__init__(atoms, indices=np.arange(atoms.get_global_number_of_atoms()), mask=None)
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.get_stress)    def get_stress(self, **kwargs: Any) -> np.ndarray:
            """Calculate stress tensor.
    
            Returns:
                np.ndarray : Stress tensor after applying filter.
            """
            stresses: np.ndarray = self.atoms.get_stress(**kwargs)
            if stresses.shape != (6,):
                raise ValueError(f"Error: stresses.shape {stresses.shape}")
            stresses[3:] = 0.0
            pressure = np.mean(stresses[:3])
            stresses[:3] = pressure
            return stresses
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.get_temperature)    def get_temperature(self) -> float:
            """Get the temperature in Kelvin.
    
            Returns:
                float : The temperature.
            """
            temperature: float = self.atoms.get_temperature()
            return temperature
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.get_total_energy)    def get_total_energy(self) -> float:
            """Get the total energy.
    
            Returns:
                float : The total energy.
            """
            total_energy: float = self.atoms.get_total_energy()
            return total_energy
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.get_kinetic_energy)    def get_kinetic_energy(self) -> float:
            """Get the kinetic energy.
    
            Returns:
                float : The kinetic energy.
            """
            kinetic_energy: float = self.atoms.get_kinetic_energy()
            return kinetic_energy
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.get_global_number_of_atoms)    def get_global_number_of_atoms(self) -> int:
            """Returns the global number of atoms in a distributed-atoms parallel simulation.
    
            Returns:
                int : The global number of atoms.
            """
            number_of_atoms: int = self.atoms.get_global_number_of_atoms()
            return number_of_atoms
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicFilter.set_cell)    def set_cell(self, cell: np.ndarray, scale_atoms: bool = False, apply_constraint: bool = True) -> None:
            """Set unit cell vectors.
    
            Args:
                cell (np.ndarray): Unit cell.
                scale_atoms (bool, optional):
                    Whether to fix atomic positions or move atoms with the unit cell. Defaults to
                    False.
                apply_constraint (bool, optional):
                    Whether to apply constraints to the given cell. Defaults to True.
            """
            self.atoms.set_cell(cell, scale_atoms=scale_atoms, apply_constraint=apply_constraint)
    
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicASEFilter)class IsotropicASEFilter(MatlantisASEFilter):
        """The filter that uses the IsotropicFilter implemmentation."""
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicASEFilter.__init__)    def __init__(self) -> None:
            """Initialize an instance."""
    
            def get_filter(atoms: Union[ASEAtoms, MatlantisAtoms]) -> Filter:
                if isinstance(atoms, MatlantisAtoms):
                    system = atoms.ase_atoms
                else:
                    system = atoms
    
                return IsotropicFilter(system)
    
            super(IsotropicASEFilter, self).__init__(get_filter=get_filter)
    
    
    
    
    
    [[docs]](../../matlantis_features.filters.html#matlantis_features.filters.IsotropicASEFilter.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the IsotropicASEFilter.
    
            Returns:
                dict[str, Any] : A dict containing a serialized IsotropicASEFilter.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
            }
            return res
    
    
    
