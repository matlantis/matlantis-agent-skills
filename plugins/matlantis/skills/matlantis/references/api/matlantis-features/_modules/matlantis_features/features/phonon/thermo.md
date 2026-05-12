# Source code for matlantis_features.features.phonon.thermo
    
    
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase import units
    from numpy import ndarray
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.force_constant import ForceConstantFeatureResult
    from matlantis_features.features.phonon.utils import PhononFrequency, get_k_mesh
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import (
        helmholtz_free_energy_converter,
        specific_heat_converter,
    )
    
    eps = 1e-5
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult)@dataclass
    class PostPhononThermochemistryFeatureResult:
        """A dataclass for result of PostPhononThermochemistryFeature."""
    
        kpts: List[int]
        num_modes: int
        num_effective_modes: int
        include_lattice_energy: bool
        temperature: Dict[str, List[float]]
        entropy: Dict[str, List[float]]
        specific_heat: Dict[str, List[float]]
        internal_energy: Dict[str, List[float]]
        helmholtz_free_energy: Dict[str, List[float]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plots the thermochemistry calculation result from phonon.
    
            Args:
                plt_name (str or None, optional): The file name to which the plot of thermochemistry
                    properties will be saved. Supported formats include  'png' 'jpg' 'jpeg' 'webp'
                    'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                    Defaults to None.
    
            Returns:
                go.Figure : The plot of thermochemical properties as a plotly.graph_objects.Figure instance.
    
            """
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=self.temperature["K"], y=self.entropy["J/K/g"], name="S (J/K/g)"))
            fig.add_trace(go.Scatter(x=self.temperature["K"], y=self.specific_heat["J/K/g"], name="Cv (J/K/g)"))
            fig.add_trace(go.Scatter(x=self.temperature["K"], y=self.internal_energy["kJ/g"], name="U (kJ/g)"))
            fig.add_trace(
                go.Scatter(
                    x=self.temperature["K"],
                    y=self.helmholtz_free_energy["kJ/g"],
                    name="F (kJ/g)",
                )
            )
            fig.update_xaxes(title="temperature [K]")
            fig.update_yaxes(title="themochemistry")
            fig.update_layout(legend=dict(x=1.0, xanchor="right", y=1.0, yanchor="bottom", orientation="h"))
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Converts the thermochemical properties to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The csv file path to which the dataframe will be
                    saved. If 'None' is provided, the dataframe will not be saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains temperature, specific_heat, entropy,
                    internal_energy and helmholtz_free_energy,
            """
            df = pd.DataFrame()
            for key in [
                "temperature",
                "specific_heat",
                "entropy",
                "internal_energy",
                "helmholtz_free_energy",
            ]:
                for unit in self.dict[key]:
                    df[f"{key} ({unit})"] = self.dict[key][unit]
    
            if csv_name:
                df.to_csv(csv_name, index=False)
    
            return df
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the PostPhononThermochemistryFeatureResult.
    
            Returns:
                dict[str, Any] : Thermochemistry properties in dictionary format.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature)class PostPhononThermochemistryFeature(FeatureBase):
        """The matlantis-feature for calculating thermochemistry properties from phonon."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature.html#matlantis_features.features.phonon.thermo.PostPhononThermochemistryFeature.__init__)    def __init__(
            self,
        ) -> None:  # TODO add 'properties' parameter, e.g. properties = ["entropy", "specific heat"]
            """Initialize an instance."""
            super(PostPhononThermochemistryFeature, self).__init__()
    
    
    
        def __call__(
            self,
            force_constant: ForceConstantFeatureResult,
            kpts: List[int],
            tmin: float = 0.0,
            tmax: float = 1000.0,
            tstep: float = 10.0,
            temperatures: Optional[List[float]] = None,
            include_lattice_energy: bool = False,
            scheme: str = "mp",
            shift: Optional[np.ndarray] = None,
            cutoff_frequency: float = 0.0,
        ) -> PostPhononThermochemistryFeatureResult:
            """Calculates the thermochemistry properties from the force constant.
    
            Args:
                force_constant (ForceConstantFeatureResult): The calculation result of
                    the ForceConstantFeature.
                kpts (list[int]): The number of mesh points along each axis.
                tmin (float, optional): The lower limit of temperate range to be calculated.
                    Defaults to 0.0.
                tmax (float, optional): The upper limit of temperate range to be calculated.
                    Defaults to 1000.0.
                tstep (float, optional): Creates a evenly spaced set of temperatures in the range of
                    [tmin, tmax) at intervals of tstep, at which the thermochemical
                    properties will be calculated. Defaults to 10.0.
                temperatures (list[float] or None, optional):
                    The temperatures at which the thermochemistry properties will be calculated. If
                    None is provided, the calculation temperatures will be generated according to
                    'tmin', 'tmax' and 'tstep'. Defaults to None.
                include_lattice_energy (bool, optional): Add the potential energy of the input
                    structure to the result of internal energy and helmholtz free energy. Defaults
                    to False.
                scheme (str, optional): The type of mesh grid: either 'mp' or 'gamma'. The mesh grid is
                    determined by the Monkhorst-Pack scheme ('mp') or Gamma-centered scheme ('gamma').
                    Defaults to 'mp'.
                shift (np.ndarray or None, optional): An array of shape (3,) containing the shift of
                    the mesh grid in the direction along the corresponding reciprocal axes. Defaults
                    to None.
                cutoff_frequency (float, optional): The vibration modes with frequency smaller than
                    'cutoff_frequency' will not be used in the thermochemical property calculation.
                    Imaginary frequencies will be regarded as negative. Defaults to 0.0.
            Returns:
                PostPhononThermochemistryFeatureResult :
                    The calculation results of thermochemical properties from phonons.
            """
            if temperatures is None:
                temperatures = np.arange(tmin, tmax + eps, tstep).tolist()
    
            pf = PhononFrequency(
                force_constant=force_constant.force_constant,
                supercell=force_constant.supercell,
                unit_cell_atoms=force_constant.unit_cell_atoms,
            )
    
            kpts_mh = get_k_mesh(kpts, scheme, shift)  # type: ignore
            frequency = pf.get_frequency(kpts_mh)
            frequency_effective = frequency[(frequency > cutoff_frequency) & (frequency != 0.0)]
    
            s: List[float] = []
            cv: List[float] = []
            u: List[float] = []
            a: List[float] = []
            masses = pf.mass
            potential_energy = force_constant.unit_cell_potential_energy
            for t in temperatures:
                if t < 0:
                    raise ValueError(f"Illegal temperature ({t} K). ")
                elif t < 5:
                    t = 5
                s.append(
                    _entropy(
                        t,
                        frequency_effective,
                        n_kpts=len(kpts_mh),
                        mass=masses,
                        potential_energy=potential_energy,
                        include_lattice_energy=include_lattice_energy,
                    )
                )
                cv.append(
                    _specific_heat(
                        t,
                        frequency_effective,
                        n_kpts=len(kpts_mh),
                        mass=masses,
                        potential_energy=potential_energy,
                        include_lattice_energy=include_lattice_energy,
                    )
                )
                u.append(
                    _internal_energy(
                        t,
                        frequency_effective,
                        n_kpts=len(kpts_mh),
                        mass=masses,
                        potential_energy=potential_energy,
                        include_lattice_energy=include_lattice_energy,
                    )
                )
                a.append(
                    _helmholtz_free_energy(
                        t,
                        frequency_effective,
                        n_kpts=len(kpts_mh),
                        mass=masses,
                        potential_energy=potential_energy,
                        include_lattice_energy=include_lattice_energy,
                    )
                )
    
            S = specific_heat_converter(s)
            Cv = specific_heat_converter(cv)
            U = helmholtz_free_energy_converter(u)
            A = helmholtz_free_energy_converter(a)
    
            result = PostPhononThermochemistryFeatureResult(
                **{
                    "kpts": kpts,
                    "num_modes": frequency.size,
                    "num_effective_modes": frequency_effective.size,
                    "include_lattice_energy": include_lattice_energy,
                    "temperature": {"K": temperatures},
                    "entropy": S,
                    "specific_heat": Cv,
                    "internal_energy": U,
                    "helmholtz_free_energy": A,
                }
            )
    
            return result
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): The number of atoms.
            Returns:
                FeatureCost :
                    A dataclass containing cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=10.0)
    
    
    
    
    def _entropy(
        t: float,
        frequency: ndarray,
        n_kpts: int,
        mass: ndarray,
        potential_energy: float,
        include_lattice_energy: bool,
    ) -> float:
        """Calculates the entropy at a given temperature.
    
        Args:
            t (float): The calculation temperature.
            frequency (ndarray): The phonon frequencies.
            n_kpts (int): Total number of k points in the mesh grid.
            mass (ndarray): The mass of each atom.
            potential_energy (float): The potential energy of input structure.
            include_lattice_energy (bool): Whether to include the potential energy in the calculation
                result.
        Returns:
            float : The entropy.
        """
        x = 1e-3 * frequency / 2 / units.kB / t
        S: float = (
            units.kB
            * np.sum(
                (1 / (np.exp(2 * x) - 1) + 1) * np.log(1 / (np.exp(2 * x) - 1) + 1)
                - (1 / (np.exp(2 * x) - 1)) * np.log(1 / (np.exp(2 * x) - 1))
            )
            / n_kpts
            / mass.sum()
        ).tolist()
        return S
    
    
    def _specific_heat(
        t: float,
        frequency: ndarray,
        n_kpts: int,
        mass: ndarray,
        potential_energy: float,
        include_lattice_energy: bool,
    ) -> float:
        """Calculate the specific heat at a given temperature.
    
        Args:
            t (float): The calculation temperature.
            frequency (ndarray): The phonon frequencies.
            n_kpts (int): Total number of k points in the mesh grid.
            mass (ndarray): The mass of each atom.
            potential_energy (float): The potential energy of the input structure.
            include_lattice_energy (bool): Whether to include the potential energy in the calculation
                result.
        Returns:
            float : The specific heat C_v.
        """
        x = 1e-3 * frequency / units.kB / t
        Cv: float = (units.kB * np.sum(x**2 * (np.exp(x) / (np.exp(x) - 1) ** 2)) / n_kpts / mass.sum()).tolist()
        return Cv
    
    
    def _internal_energy(
        t: float,
        frequency: ndarray,
        n_kpts: int,
        mass: ndarray,
        potential_energy: float,
        include_lattice_energy: bool,
    ) -> float:
        """Calculates the internal energy at a given temperature.
    
        Args:
            t (float): The calculation temperature.
            frequency (ndarray): The phonon frequencies.
            n_kpts (int): Total number of k points in the mesh grid.
            mass (ndarray): The mass of each atom.
            potential_energy (float): The potential energy of the input structure.
            include_lattice_energy (bool): Whether to include the potential energy in the calculation
                result.
        Returns:
            float : The internal energy.
        """
        U: float = (
            np.sum((1e-3 * frequency / 2) + 1e-3 * frequency / (np.exp(1e-3 * frequency / units.kB / t) - 1)) / n_kpts / mass.sum()
        ).tolist()
    
        if include_lattice_energy:
            return float(U + potential_energy / mass.sum())  # unit eV/amu
        else:
            return U
    
    
    def _helmholtz_free_energy(
        t: float,
        frequency: ndarray,
        n_kpts: int,
        mass: ndarray,
        potential_energy: float,
        include_lattice_energy: bool,
    ) -> float:
        """Calculates the Helmholtz free energy at a given temperature.
    
        Args:
            t (float): The calculation temperature.
            frequency (ndarray): The phonon frequencies.
            n_kpts (int): Total number of k points in the mesh grid.
            mass (ndarray): The mass of each atom.
            potential_energy (float): The potential energy of the input structure.
            include_lattice_energy (bool): Whether to include the potential energy in the calculation
                result.
        Returns:
            float : The Helmholtz free energy.
        """
        x = 1e-3 * frequency / units.kB / t
        A: float = (units.kB * t * np.sum((x / 2) + np.log(1 - np.exp(-x))) / n_kpts / mass.sum()).tolist()
        if include_lattice_energy:
            return float(A + potential_energy / mass.sum())  # unit eV/amu
        else:
            return A
    
