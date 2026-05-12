# Source code for matlantis_features.features.common.gas_formation_enthalpy
    
    
    import json
    import os
    import warnings
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional, Union
    
    import numpy as np
    import pandas as pd
    from ase import Atoms, units
    from ase.build import molecule
    from ase.cell import Cell
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.common.gas_thermo import PostVibrationGasThermoFeature
    from matlantis_features.features.common.opt import FireASEOptFeature, OptFeatureBase
    from matlantis_features.features.common.vibration import VibrationFeature
    from matlantis_features.features.phonon.force_constant import ForceConstantFeature
    from matlantis_features.features.phonon.thermo import PostPhononThermochemistryFeature
    from matlantis_features.utils.calculators import get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    from matlantis_features.utils.units_util import helmholtz_free_energy_converter, pressure_converter
    
    with open(os.path.dirname(__file__) + "/standard_formation_enthalpy_parameters.json", "r") as f:
        standard_enthalpy_parameters = json.load(f)
    
    calc_mode_alias_mapping = {
        "crystal_u0": "pbe",
        "crystal_u0_plus_d3": "pbe_plus_d3",
        "crystal": "pbe_u",
        "crystal_plus_d3": "pbe_u_plus_d3",
        "molecule": "wb97xd",
    }
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult)@dataclass
    class ComplexGasFormationEnthalpyFeatureResult:
        """A dataclass for result of ComplexFormationEnthalpyFeature."""
    
        temperature: Dict[str, List[float]]
        pressure: Dict[str, List[float]]
        formation_enthalpy: Dict[str, List[float]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Converts the formation enthalpy calculation result to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The csv file path to which the dataframe will be
                    saved. If 'None' is provided, the dataframe will not be saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains temperature, pressure and formation_enthalpy.
            """
            df = pd.DataFrame()
            for key in [
                "temperature",
                "pressure",
                "formation_enthalpy",
            ]:
                for unit in self.dict[key]:
                    df[f"{key} ({unit})"] = self.dict[key][unit]
    
            if csv_name:
                df.to_csv(csv_name, index=False)
    
            return df
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the ComplexFormationEnthalpyFeatureResult.
    
            Returns:
                dict[str, Any] : Formation enthalpy in dictionary format.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature)class ComplexGasFormationEnthalpyFeature(FeatureBase):
        """The matlantis-feature for calculating gas formation enthalpy."""
    
        allowed_elements = np.array([1, 6, 7, 8])
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature.html#matlantis_features.features.common.gas_formation_enthalpy.ComplexGasFormationEnthalpyFeature.__init__)    def __init__(self, opt_feature: OptFeatureBase, vibration_feature: VibrationFeature) -> None:
            """Initialize an instance.
    
            .. note::
                Since ``ComplexGasFormationEnthalpyFeature`` instantiates an estimator in `__init__`,
                it won't work with multiprocessing parallelization. If you want to use this feature
                in multiple processes, please make sure you use `loky` for the backend.
                For more information, please refer to the `Matlantis guidebook </api/resource/documents/
                matlantis-guidebook/en/notebooks.html#the-use-of-estimator-in-multi-process-computations>`_.
    
            Args:
                opt_feature (OptFeatureBase):
                  Optimizer feature object used for structure optimization before the enthalpy
                  calculation.
                vibration_feature (VibrationFeature):
                  Feature object for calculating molecule vibration.
            """
            super(ComplexGasFormationEnthalpyFeature, self).__init__()
            with self.init_scope():
                self.opt_feature = opt_feature
                self.vib_feature = vibration_feature
                self.post_vib_feature = PostVibrationGasThermoFeature()
                common_estimator_fn = self._validate_estimator_fn(opt_feature, vibration_feature)
                self._opt_feature_solid = FireASEOptFeature(n_run=1000, filter=True, estimator_fn=common_estimator_fn)
                self._fc_feature = ForceConstantFeature(supercell=(3, 3, 3), delta=0.1, estimator_fn=common_estimator_fn)
                self._post_phonon_feature = PostPhononThermochemistryFeature()
    
    
    
        def __call__(
            self,
            atoms: Union[ASEAtoms, MatlantisAtoms],
            temperatures: List[float],
            pressures: List[float],
        ) -> ComplexGasFormationEnthalpyFeatureResult:
            """Calculate the formation enthalpy of the gas molecule.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms): The input structure.
                temperatures (list[float]):
                    The temperatures at which the formation enthalpy will be calculated.
                pressures (list[float]):
                    The pressures at which the formation enthalpy will be calculated. The unit of
                    pressure is eV/Angstrom^3.
            Returns:
                ComplexGasFormationEnthalpyFeatureResult :
                    The calculation results of formation enthalpy.
            """
            if isinstance(atoms, MatlantisAtoms):
                system = atoms.ase_atoms
            else:
                system = atoms
    
            assert np.all(np.isin(atoms.numbers, self.allowed_elements)), "Only C H O N elements are supported."
            enthalpy = self._get_gas_enthalpy(atoms, temperatures, pressures)
    
            if np.any(system.get_pbc()):
                warnings.warn("The input structure is supposed to be a gas molecule without pbc.", UserWarning)
    
            ref_enthalpy = np.zeros(len(temperatures))
            if 1 in system.numbers:
                n_H = np.sum(system.numbers == 1)
                ref_enthalpy += self._get_H_enthalpy(temperatures, pressures) * n_H
            if 6 in system.numbers:
                n_C = np.sum(system.numbers == 6)
                ref_enthalpy += self._get_C_enthalpy(temperatures, pressures) * n_C
            if 7 in system.numbers:
                n_N = np.sum(system.numbers == 7)
                ref_enthalpy += self._get_N_enthalpy(temperatures, pressures) * n_N
            if 8 in system.numbers:
                n_O = np.sum(system.numbers == 8)
                ref_enthalpy += self._get_O_enthalpy(temperatures, pressures) * n_O
    
            formation_enthalpy = (enthalpy - ref_enthalpy).tolist()
    
            return ComplexGasFormationEnthalpyFeatureResult(
                temperature={"K": temperatures},
                pressure=pressure_converter(pressures),
                formation_enthalpy=helmholtz_free_energy_converter(formation_enthalpy, per_mass=False),
            )
    
        def _get_gas_enthalpy(self, atoms: ASEAtoms, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            opt_result = self.opt_feature(atoms)
            vib_result = self.vib_feature(opt_result.atoms)
            thermo_result = self.post_vib_feature(vib_result, temperatures, pressures)
            return np.array(thermo_result.enthalpy["eV"])  # type: ignore
    
        def _get_H_enthalpy(self, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            atoms = molecule("H2")
            return self._get_gas_enthalpy(atoms, temperatures, pressures) / 2.0  # type: ignore
    
        def _get_O_enthalpy(self, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            atoms = molecule("O2")
            return self._get_gas_enthalpy(atoms, temperatures, pressures) / 2.0  # type: ignore
    
        def _get_N_enthalpy(self, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            atoms = molecule("N2")
            return self._get_gas_enthalpy(atoms, temperatures, pressures) / 2.0  # type: ignore
    
        def _get_solid_enthalpy(self, atoms: Atoms, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            volume = atoms.get_volume()
            mass = atoms.get_masses().sum()
    
            opt_result = self._opt_feature_solid(atoms)
            fc_result = self._fc_feature(opt_result.atoms)
            thermo_result = self._post_phonon_feature(fc_result, kpts=[10, 10, 10], temperatures=temperatures, include_lattice_energy=True)
    
            internal_energy: np.ndarray = np.array(thermo_result.internal_energy["eV/amu"]) * mass
            return internal_energy + volume * np.array(pressures)  # type: ignore
    
        def _get_C_enthalpy(self, temperatures: List[float], pressures: List[float]) -> np.ndarray:
            graphite = Atoms(
                scaled_positions=[
                    [0, 0, 0.25],
                    [0, 0, 0.75],
                    [0.33333, 0.66667, 0.25],
                    [0.66667, 0.33333, 0.75],
                ],
                cell=Cell.fromcellpar([2.464, 2.464, 6.711, 90.0, 90.0, 120.0]),
                numbers=[6, 6, 6, 6],
                pbc=True,
            )
            return self._get_solid_enthalpy(graphite, temperatures, pressures) / 4  # type: ignore
    
        def _validate_estimator_fn(self, opt_feature: OptFeatureBase, vib_feature: VibrationFeature) -> Optional[EstimatorFnType]:
            """
            Validate that the estimator functions for the optimization and vibration features are consistent.
    
            Args:
                opt_feature (OptFeatureBase): Optimization feature object.
                vib_feature (VibrationFeature): Vibration feature object.
    
            Returns:
                Optional[EstimatorFnType]: The common estimator function if both features have the same estimator_fn,
                otherwise None.
    
            Raises:
                RuntimeError: If only one of the features has an estimator function set
                or if the estimator functions are different.
            """
            if opt_feature.estimator_fn is not None and vib_feature.estimator_fn is not None:
                if opt_feature.estimator_fn != vib_feature.estimator_fn:
                    raise RuntimeError(
                        "Different `estimator_fn` are passed to vibration feature and opt feature. "
                        "You must pass the same `estimator_fn` to use the same estimator across "
                        "vibration feature and opt feature."
                    )
                else:
                    return opt_feature.estimator_fn
            elif (opt_feature.estimator_fn is None and vib_feature.estimator_fn is not None) or (
                opt_feature.estimator_fn is not None and vib_feature.estimator_fn is None
            ):
                raise RuntimeError(
                    "If `estimator_fn` is set in either vibration feature or opt feature, the same estimator must be set in "
                    "the other feature as well."
                )
            else:
                return None
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeatureResult.html#matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeatureResult)@dataclass
    class GasStandardFormationEnthalpyFeatureResult:
        """A dataclass for result of GasStandardFormationEnthalpyFeature."""
    
        standard_formation_enthalpy: Dict[str, List[float]]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature.html#matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature)class GasStandardFormationEnthalpyFeature(FeatureBase):
        """The matlantis-feature for calculating standard formation enthalpy of gas molecule."""
    
        parameters = standard_enthalpy_parameters
        allowed_elements = np.array([1, 6, 7, 8])
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature.html#matlantis_features.features.common.gas_formation_enthalpy.GasStandardFormationEnthalpyFeature.__init__)    def __init__(self, opt_feature: OptFeatureBase, vibration_feature: VibrationFeature) -> None:
            """Initialize an instance.
    
            .. note::
                Since ``GasStandardFormationEnthalpyFeature`` instantiates an estimator in `__init__`,
                it won't work with multiprocessing parallelization. If you want to use this feature
                in multiple processes, please make sure you use `loky` for the backend.
                For more information, please refer to the `Matlantis guidebook </api/resource/documents/
                matlantis-guidebook/en/notebooks.html#the-use-of-estimator-in-multi-process-computations>`_.
    
            Args:
                opt_feature (OptFeatureBase):
                  Optimizer feature object used for structure optimization before the enthalpy
                  calculation.
                vibration_feature (VibrationFeature):
                  Feature object for calculating molecule vibration.
            """
            super(GasStandardFormationEnthalpyFeature, self).__init__()
            with self.init_scope():
                self.opt_feature = opt_feature
                self.vib_feature = vibration_feature
                self.post_vib_feature = PostVibrationGasThermoFeature()
    
                # NOTE (himkt)
                # ``estimator_fn`` hides ``model_version`` and ``calc_mode``.
                # Here it instantiates estimators and check if both have the
                # same ``model_version`` and ``calc_mode``.
                if self.vib_feature.estimator_fn is not None:
                    if self.opt_feature.estimator_fn is None:
                        raise RuntimeError(
                            "If `estimator_fn` is set in vibration feature, the same estimator must be set in opt feature as well."
                        )
    
                    vib_estimator = self.vib_feature.estimator_fn()
                    opt_estimator = self.opt_feature.estimator_fn()
    
                    if not isinstance(vib_estimator, Estimator):
                        raise RuntimeError("Currently, PFP estimator is only supported.")
    
                    if vib_estimator.model_version != opt_estimator.model_version or vib_estimator.calc_mode != opt_estimator.calc_mode:
                        raise RuntimeError(
                            "Different ``estimator_fn`` are passed to"
                            "vibration feature and opt feature. "
                            "You must pass the same ``estimator_fn`` to "
                            "use the same estimator across vibration feature "
                            "and opt feature."
                        )
    
                    model_version = vib_estimator.model_version
                    calc_mode = vib_estimator.calc_mode
                    if isinstance(calc_mode, EstimatorCalcMode):
                        calc_mode = calc_mode.value
                else:
                    # model_version = os.environ.get("MATLANTIS_PFP_MODEL_VERSION", "latest").lower()
                    # calc_mode = os.environ.get("MATLANTIS_PFP_CALC_MODE", "CRYSTAL").lower()
                    estimator = get_calculator().estimator
                    model_version = estimator.model_version
                    if estimator.calc_mode is None:
                        calc_mode = "crystal"  # TODO: fix this hard-coded value
                    else:
                        calc_mode = estimator.calc_mode.name.lower()
    
                self.model_version = model_version
                self.calc_mode = calc_mode
    
    
    
        def __call__(
            self,
            atoms: Union[ASEAtoms, MatlantisAtoms],
        ) -> GasStandardFormationEnthalpyFeatureResult:
            """Calculate the formation enthalpy of the gas molecule.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms): The input structure.
            Returns:
                GasStandardFormationEnthalpyFeatureResult :
                    The calculation results of standard formation enthalpy.
            """
            if isinstance(atoms, MatlantisAtoms):
                system = atoms.ase_atoms
            else:
                system = atoms
    
            assert np.all(np.isin(atoms.numbers, self.allowed_elements)), "Only C H O N elements are supported."
    
            if self.model_version not in self.parameters:
                raise ValueError(
                    f"The model version '{self.model_version}' is not supported in the "
                    "GasStandardFormationEnthalpyFeature. Supported model versions are "
                    + " ".join([f"'{k}'" for k in self.parameters.keys()])
                    + "."
                )
            calc_mode_ = calc_mode_alias_mapping.get(self.calc_mode, self.calc_mode)
            if calc_mode_ not in self.parameters[self.model_version]:
                raise ValueError(
                    f"The calculation mode '{calc_mode_}' is not supported in the "
                    "GasStandardFormationEnthalpyFeature. Supported calculation modes are "
                    + " ".join([f"'{k}'" for k in self.parameters[self.model_version].keys()])
                    + "."
                )
    
            enthalpy = self._get_gas_enthalpy(atoms)
    
            scaler = 1000 * units.J / units._Nav
            ref_enthalpy = np.array([0.0])
            if 1 in system.numbers:
                n_H = np.sum(system.numbers == 1)
                ref_enthalpy += self.parameters[self.model_version][calc_mode_]["H_H"] * n_H
            if 6 in system.numbers:
                n_C = np.sum(system.numbers == 6)
                ref_enthalpy += self.parameters[self.model_version][calc_mode_]["H_C"] * n_C
            if 7 in system.numbers:
                n_N = np.sum(system.numbers == 7)
                ref_enthalpy += self.parameters[self.model_version][calc_mode_]["H_N"] * n_N
            if 8 in system.numbers:
                n_O = np.sum(system.numbers == 8)
                ref_enthalpy += self.parameters[self.model_version][calc_mode_]["H_O"] * n_O
    
            formation_enthalpy = (enthalpy - ref_enthalpy * scaler).tolist()
    
            return GasStandardFormationEnthalpyFeatureResult(
                standard_formation_enthalpy=helmholtz_free_energy_converter(formation_enthalpy, per_mass=False),
            )
    
        def _get_gas_enthalpy(self, atoms: ASEAtoms) -> np.ndarray:
            opt_result = self.opt_feature(atoms)
            vib_result = self.vib_feature(opt_result.atoms)
            thermo_result = self.post_vib_feature(vib_result, temperatures=[298.15], pressures=[6.32421e-07])
            return np.array(thermo_result.enthalpy["eV"])  # type: ignore
    
    
    
