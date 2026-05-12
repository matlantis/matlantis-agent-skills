# Source code for matlantis_features.features.common.gas_thermo
    
    
    import warnings
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional
    
    import pandas as pd
    from ase import units
    from ase.thermochemistry import IdealGasThermo
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.common.vibration import VibrationFeatureResult
    from matlantis_features.utils.atoms_util import (
        MoleculeLinearity,
        check_linearity,
        check_symmetry_number,
    )
    from matlantis_features.utils.units_util import (
        helmholtz_free_energy_converter,
        pressure_converter,
        specific_heat_converter,
    )
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult.html#matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult)@dataclass
    class PostVibrationGasThermoFeatureResult:
        """A dataclass for result of PostVibrationGasThermoFeature."""
    
        temperature: Dict[str, List[float]]
        pressure: Dict[str, List[float]]
        zero_point_energy: float
        enthalpy: Dict[str, List[float]]
        entropy: Dict[str, List[float]]
        gibbs_free_energy: Dict[str, List[float]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult.html#matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Converts the thermochemical properties to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The csv file path to which the dataframe will be
                    saved. If 'None' is provided, the dataframe will not be saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains temperature, pressure, entropy, enthalpy
                    and helmholtz_free_energy,
            """
            df = pd.DataFrame()
            for key in [
                "temperature",
                "pressure",
                "entropy",
                "enthalpy",
                "gibbs_free_energy",
            ]:
                for unit in self.dict[key]:
                    df[f"{key} ({unit})"] = self.dict[key][unit]
    
            if csv_name:
                df.to_csv(csv_name, index=False)
    
            return df
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the PostVibrationGasThermoFeatureResult.
    
            Returns:
                dict[str, Any] : Thermochemistry properties in dictionary format.
            """
            return asdict(self)
    
    
    
    
    def _check_vibration_frequency(vibration_results: VibrationFeatureResult) -> bool:
        imaginary = False
        ase_atoms = vibration_results.atoms.ase_atoms
        geometry = check_linearity(ase_atoms)
        frequencies = vibration_results.frequencies
        image_frequencies = vibration_results.image_frequencies
    
        if geometry is MoleculeLinearity.nonlinear:
            invalid_modes = 6
        elif geometry is MoleculeLinearity.linear:
            invalid_modes = 5
    
        for i in range(invalid_modes, len(frequencies)):
            if frequencies[i] <= 0.0:
                warnings.warn(
                    f"Imaginary frequency is found at the {i} th vibration mode of a {geometry.name} "
                    f"molecule. The imaginary frequency is {frequencies[i]} + {image_frequencies[i]}i"
                )
                imaginary = True
        return imaginary
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeature.html#matlantis_features.features.common.gas_thermo.PostVibrationGasThermoFeature)class PostVibrationGasThermoFeature(FeatureBase):
        """The matlantis-feature for calculating the thermal property of gas molecule from the \
            molecular vibration. The ideal-gas-approximation is used.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/vibration_thermo.html>`_.
        """
    
        def __call__(
            self,
            vibration_results: VibrationFeatureResult,
            temperatures: List[float],
            pressures: List[float],
        ) -> PostVibrationGasThermoFeatureResult:
            """Initialize an instance.
    
            Args:
                vibration_results (VibrationFeatureResult): The calculation result of VibrationFeature.
                temperatures (list[float]):
                    The temperatures at which the thermochemistry properties will be calculated.
                pressures (list[float]):
                    The pressures at which the thermochemistry properties will be calculated. The unit
                    of pressure is eV/Angstrom^3.
            Returns:
                PostVibrationGasThermoFeatureResult :
                The calculation results including zero point energy of molecule, enthalpy, entropy
                gibbs free energy at the given temperatures and pressures.
            """
            ase_atoms = vibration_results.atoms.ase_atoms
            geometry = check_linearity(ase_atoms)
            symmetry_number = check_symmetry_number(ase_atoms)
    
            imaginary = _check_vibration_frequency(vibration_results)
    
            frequencies = vibration_results.frequencies.copy()
            if imaginary:
                # Imaginary frequency was found.
                # In order to prevent the error caused by imaginary frequencies, we add a small number.
                frequencies[frequencies == 0.0] += 1e-4
                warnings.warn(
                    "Imaginary frequencies are found. Please check the vibration calculation results. "
                    "The calculation results of thermal properties, especially entropy(S) and Gibbs "
                    "free energy(G), are ill defined due to the imagenary frequency."
                )
    
            ideal_gas_thermo = IdealGasThermo(
                vib_energies=frequencies,
                geometry=geometry.name,
                potentialenergy=vibration_results.potential_energy,
                atoms=ase_atoms,
                symmetrynumber=symmetry_number,
                spin=0,
            )
    
            zpe = ideal_gas_thermo.get_ZPE_correction()
            enthalpy = []
            entropy = []
            gibbs = []
            for t, p in zip(temperatures, pressures):
                enthalpy.append(ideal_gas_thermo.get_enthalpy(t, verbose=False))
                entropy.append(ideal_gas_thermo.get_entropy(t, p / units.Pascal, verbose=False))
                gibbs.append(ideal_gas_thermo.get_gibbs_energy(t, p / units.Pascal, verbose=False))
    
            return PostVibrationGasThermoFeatureResult(
                temperature={"K": temperatures},
                pressure=pressure_converter(pressures),
                zero_point_energy=zpe,
                enthalpy=helmholtz_free_energy_converter(enthalpy, per_mass=False),
                entropy=specific_heat_converter(entropy, per_mass=False),
                gibbs_free_energy=helmholtz_free_energy_converter(gibbs, per_mass=False),
            )
    
    
    
