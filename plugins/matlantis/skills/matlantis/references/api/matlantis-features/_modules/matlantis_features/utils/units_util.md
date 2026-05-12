# Source code for matlantis_features.utils.units_util
    
    
    from typing import Dict, List, Union
    
    import numpy as np
    from ase import units
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.specific_heat_converter.html#matlantis_features.utils.units_util.specific_heat_converter)def specific_heat_converter(specific_heat: List[float], per_mass: bool = True) -> Dict[str, List[float]]:
        """Change the unit of specific heat (or entropy) from eV/K or eV/K/amu.
    
        Args:
            specific_heat (List[float]): unit eV/K/amu
    
        Returns:
            Dict[str, List[float]]: the specific heat (or entropy) to 'eV/K/amu', 'J/K/g' and 'kJ/K/g'.
        """
        if per_mass:
            return {
                "eV/K/amu": specific_heat,
                "J/K/g": [i * units.kg / units.J / 1000 for i in specific_heat],
                "kJ/K/g": [i * units.kg / units.J / 1e6 for i in specific_heat],
                "kcal/K/g": [i * units.kg / units.kcal / 1000 for i in specific_heat],
            }
        else:
            return {
                "eV/K": specific_heat,
                "J/K/mol": [i / units.J * units._Nav for i in specific_heat],
                "kJ/K/mol": [i / units.J * units._Nav / 1000 for i in specific_heat],
                "kcal/K/mol": [i / units.kcal * units._Nav for i in specific_heat],
            }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.thermal_expansion_converter.html#matlantis_features.utils.units_util.thermal_expansion_converter)def thermal_expansion_converter(thermal_expansion: List[float]) -> Dict[str, List[float]]:
        """add the units for volumetric expansion coefficient and linear expansion coefficient.
    
        Args:
            thermal_expansion (List[float]): [description]
    
        Returns:
            Dict[str, List[float]]: the thermal expansion coefficient in '1/K'.
        """
        return {
            "1/K": thermal_expansion,
        }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.helmholtz_free_energy_converter.html#matlantis_features.utils.units_util.helmholtz_free_energy_converter)def helmholtz_free_energy_converter(helmholtz_free_energy: List[float], per_mass: bool = True) -> Dict[str, List[float]]:
        """Change the unit of helmholtz free energy (or internal energy) from eV or eV/amu.
    
        Args:
            helmholtz_free_energy (List[float]): unit eV/amu
    
        Returns:
            Dict[str, List[float]]: the specific heat (or entropy) to 'eV/amu', 'J/g' and 'kJ/g'.
        """
        if per_mass:
            return {
                "eV/amu": helmholtz_free_energy,
                "J/g": [i * units.kg / units.J / 1000 for i in helmholtz_free_energy],
                "kJ/g": [i * units.kg / units.J / 1e6 for i in helmholtz_free_energy],
                "kcal/g": [i * units.kg / units.kcal / 1000 for i in helmholtz_free_energy],
            }
        else:
            return {
                "eV": helmholtz_free_energy,
                "J/mol": [i / units.J * units._Nav for i in helmholtz_free_energy],
                "kJ/mol": [i / units.J * units._Nav / 1e3 for i in helmholtz_free_energy],
                "kcal/mol": [i / units.kcal * units._Nav for i in helmholtz_free_energy],
            }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.viscosity_converter.html#matlantis_features.utils.units_util.viscosity_converter)def viscosity_converter(viscosity: List[float]) -> Dict[str, List[float]]:
        """Convert viscosity values from amu/Angstrom/fs to SI units.
    
        Args:
            viscosity (List[float]): Viscosity values in amu/Angstrom/fs.
    
        Returns:
            Dict[str, List[float]]:
        """
        return {
            "amu/A/fs": viscosity,
            "Pa s": [i * 1e25 / units.kg for i in viscosity],
            "mPa s": [i * 1e28 / units.kg for i in viscosity],
        }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.diffusion_coefficient_converter.html#matlantis_features.utils.units_util.diffusion_coefficient_converter)def diffusion_coefficient_converter(diffusion_coefficient: Union[float, List[float]]) -> Union[Dict[str, float], Dict[str, List[float]]]:
        if isinstance(diffusion_coefficient, List):
            return {
                "A^2/fs": [i for i in diffusion_coefficient],
                "cm^2/s": [i * 0.1 for i in diffusion_coefficient],
            }
        else:
            return {
                "A^2/fs": diffusion_coefficient,
                "cm^2/s": diffusion_coefficient * 0.1,
            }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.diffusion_coefficient_converter_old.html#matlantis_features.utils.units_util.diffusion_coefficient_converter_old)def diffusion_coefficient_converter_old(
        diffusion_coefficient: List[float],
    ) -> Dict[str, List[float]]:
        return {
            "A^2/fs": [i * units.fs for i in diffusion_coefficient],
            "cm^2/s": [i * units.fs * 0.1 for i in diffusion_coefficient],
        }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.thermal_conductivity_converter.html#matlantis_features.utils.units_util.thermal_conductivity_converter)def thermal_conductivity_converter(thermal_conductivity: List[float]) -> Dict[str, List[float]]:
        return {
            "eV/fs/Angstrom/K": thermal_conductivity,
            "J/s/m/K": [i * 1e25 / units.J for i in thermal_conductivity],
            "W/m/K": [i * 1e25 / units.J for i in thermal_conductivity],
        }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.frequency_converter.html#matlantis_features.utils.units_util.frequency_converter)def frequency_converter(frequency: List[float]) -> Dict[str, List[float]]:
        mev_to_THz = 0.001 / (units._hbar * 2 * np.pi * 1e12 * units.J)
        mev_to_cm = mev_to_THz * 1e10 / units._c
        return {
            "meV": frequency,
            "eV": [i / 1000.0 for i in frequency],
            "THz": [i * mev_to_THz for i in frequency],
            "cm^-1": [i * mev_to_cm for i in frequency],
        }
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.utils.units_util.pressure_converter.html#matlantis_features.utils.units_util.pressure_converter)def pressure_converter(pressure: List[float]) -> Dict[str, List[float]]:
        return {
            "eV/Angstrom^3": pressure,
            "Pa": [p / units.Pascal for p in pressure],
            "GPa": [p / units.GPa for p in pressure],
            "bar": [p / units.bar for p in pressure],
        }
    
    
    
