# Source code for matlantis_features.features.phonon.qha_gibbs_free_energy
    
    
    import logging
    from dataclasses import dataclass
    from typing import List, Optional, cast
    
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np
    import plotly.graph_objects as go
    from ase import units
    from ase.cell import Cell
    from plotly.subplots import make_subplots
    from scipy.interpolate import CloughTocher2DInterpolator
    
    from matlantis_features.atoms import MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.force_constant import ForceConstantFeatureResult
    from matlantis_features.features.phonon.thermo import (
        PostPhononThermochemistryFeature,
        PostPhononThermochemistryFeatureResult,
    )
    from matlantis_features.functions.eos import BirchMurnaghanEOS, EOSFitError, EOSVolumeError
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult)@dataclass
    class PostPhononQHAGibbsFeatureResult:
        """A dataclass for result of PostPhononQHAGibbsFeatureResult."""
    
        force_constant_list: List[ForceConstantFeatureResult]
        unit_cell_atoms: List[MatlantisAtoms]
        kpts: List[int]
        temperatures: List[float]
        pressures: List[float]
        thermochemistry: List[PostPhononThermochemistryFeatureResult]
        equilibrium_volume: np.ndarray
        helmholtz_free_energy: np.ndarray
        gibbs_free_energy: np.ndarray
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plots all calculation results from QHA.
    
            Args:
                plt_name (str or None, optional): The file name to which the plot of Gibbs free
                    energy will be saved. Supported formats include  'png' 'jpg' 'jpeg' 'webp'
                    'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                    Defaults to None.
    
            Returns:
                go.Figure : The plot of Gibbs free energy as a plotly.graph_objects.Figure instance.
            """
            F = np.empty([len(self.temperatures), len(self.force_constant_list)])
            V = np.empty(len(self.force_constant_list))
    
            for i, (thermochemistry, unit_cell_atoms) in enumerate(zip(self.thermochemistry, self.unit_cell_atoms)):
                F[:, i] = (
                    np.array(thermochemistry.helmholtz_free_energy["eV/amu"])
                    * unit_cell_atoms.ase_atoms.get_masses().sum()
                    / len(unit_cell_atoms.ase_atoms)
                )
                V[i] = unit_cell_atoms.ase_atoms.get_volume() / len(unit_cell_atoms.ase_atoms)
    
            if len(self.temperatures) > 5:
                n_stride_t = len(self.temperatures) // 5
            else:
                n_stride_t = 1
            range_t = np.arange(0, len(F), n_stride_t)
            if len(F) - 1 not in range_t:
                range_t = np.append(range_t, len(F) - 1)
    
            fig = make_subplots(rows=1, cols=2)
            for i in range_t:
                fig.add_trace(go.Scatter(x=V, y=F[i], text=f"{self.temperatures[i]:.1f} K", showlegend=False))
    
            if len(self.pressures) > 5:
                n_stride_p = len(self.pressures) // 5
            else:
                n_stride_p = 1
            range_p = np.arange(0, len(self.pressures), n_stride_p)
            if len(self.pressures) - 1 not in range_p:
                range_p = np.append(range_p, len(self.pressures) - 1)
    
            for i in range_p:
                fig.add_trace(
                    go.Scatter(
                        x=self.equilibrium_volume[:, i],
                        y=self.helmholtz_free_energy[:, i],
                        text=f"{self.pressures[i] / units.GPa:.1f} GPa",
                        showlegend=False,
                    )
                )
    
            for i in range_p:
                fig.add_trace(
                    go.Scatter(
                        x=self.temperatures,
                        y=self.gibbs_free_energy[:, i],
                        text=f"{self.pressures[i] / units.GPa:.1f} GPa",
                        name=f"{self.pressures[i] / units.GPa:.1f} GPa",
                        showlegend=True,
                    ),
                    row=1,
                    col=2,
                )
            fig.update_xaxes(title="Volume (A^3/atom)", row=1, col=1)
            fig.update_yaxes(title="Helmoholtz free energy (eV/atom)", row=1, col=1)
            fig.update_xaxes(title="temperature (K)", row=1, col=2)
            fig.update_yaxes(title="Gibbs free energy (eV/atom)", row=1, col=2)
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeatureResult.plot_3d)    def plot_3d(self) -> go.Figure:
            """Plots the surface of Gibbs free energy G(T, P).
    
            Returns:
                go.Figure : The plot of Gibbs free energy as a plotly.graph_objects.Figure instance.
            """
            fig = go.Figure()
            fig.add_trace(
                go.Surface(
                    z=self.gibbs_free_energy,
                    y=np.array(self.temperatures),
                    x=np.array(self.pressures) / units.GPa,
                )
            )
            fig.update_layout(
                scene=dict(xaxis_title="P (GPa)", yaxis_title="T (K)", zaxis_title="G (eV/atom)"),
            )
            return fig
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature)class PostPhononQHAGibbsFeature(FeatureBase):
        """The matlantis-feature for calculating Gibbs free energy from Quasi-harmonic approximation."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature.html#matlantis_features.features.phonon.qha_gibbs_free_energy.PostPhononQHAGibbsFeature.__init__)    def __init__(self) -> None:
            """Initialize an instance."""
            super(PostPhononQHAGibbsFeature, self).__init__()
            self.thermochemistry = PostPhononThermochemistryFeature()
    
    
    
        def __call__(
            self,
            force_constant_list: List[ForceConstantFeatureResult],
            kpts: List[int],
            temperatures: List[float],
            pressures: List[float],
        ) -> PostPhononQHAGibbsFeatureResult:
            """Calculates the Gibbs free energy from a list of force constants.
    
            Args:
                force_constant_list (list[ForceConstantFeatureResult]): A list of force constant.
                    The force constants are obtained by the ForceConstantFeature.
                    The unit cell should be same except the volume is different.
                    The `len(force_constant_list)` must be larger than 3.
                kpts  (list[int]):
                    The number of mesh points along each axis
                    in the Helmholtz free energy calculation
                temperatures (list[float]):
                    The temperatures at which the Gibbs free energy will be calculated.
                pressures (list[float]):
                    The pressures at which the Gibbs free energy will be calculated.
            Returns:
                PostPhononQHAGibbsFeatureResult :
                    The calculation results of Gibbs free energy calculation.
            """
            logger = logging.getLogger(__name__)
            if len(force_constant_list) < 4:
                raise ValueError(
                    "Please input more than three 'force_constant' into the 'force_constant_list', "
                    "The equilibrium volume at specific temperatures is calculated from the EOS, "
                    "but 'RankingWarning' will occur during the EOS fitting "
                    "if the number of 'force_constant' is too small."
                )
            # calculation of helmholtz free energy
            a = []
            v = []
            unit_cell_atoms = []
            thermochemistry = []
    
            for force_constant in force_constant_list:
                thermo_results: PostPhononThermochemistryFeatureResult = self.thermochemistry(
                    force_constant, kpts, temperatures=temperatures, include_lattice_energy=True
                )
                atoms = force_constant.unit_cell_atoms.ase_atoms
                a.append(np.array(thermo_results.helmholtz_free_energy["eV/amu"]) * atoms.get_masses().sum() / len(atoms))
                v.append(Cell(atoms.cell).volume / len(atoms))
                unit_cell_atoms.append(force_constant.unit_cell_atoms)
                thermochemistry.append(thermo_results)
    
            A = np.array(a)
            V = np.array(v)
    
            # fitting of eos
            equilibrium_volume = np.full([len(temperatures), len(pressures)], np.nan)
            helmholtz_free_energy = np.full([len(temperatures), len(pressures)], np.nan)
            gibbs_free_energy = np.full([len(temperatures), len(pressures)], np.nan)
    
            for i, t in enumerate(temperatures):
                try:
                    eos = BirchMurnaghanEOS.from_fit(V, A[:, i])
                except EOSFitError:
                    logger.warn("Cannot fit the EOS at %f K" % t)
                else:
                    for j, p in enumerate(pressures):
                        try:
                            V_p = eos.get_volume_under_pressure(p)
                        except EOSVolumeError:
                            logger.warn(f"Cannot get equlibrium volume at {t} K {p / units.GPa} GPa")
                        else:
                            if V_p == 0.0 or V_p > max(v) or V_p < min(v):
                                logger.warn(
                                    f"The equlibrium volume {V_p:.1f} exceed valid range "
                                    f"({min(v):.1f} {max(v):.1f}) at {t:.2f} K {p / units.GPa:.2f} GPa"
                                )
                            else:
                                equilibrium_volume[i, j] = V_p
                                helmholtz_free_energy[i, j] = cast(np.ndarray, eos(V_p))
                                gibbs_free_energy[i, j] = cast(np.ndarray, eos(V_p) + p * V_p)
    
            result = PostPhononQHAGibbsFeatureResult(
                force_constant_list=force_constant_list,
                unit_cell_atoms=unit_cell_atoms,
                kpts=kpts,
                temperatures=temperatures,
                pressures=pressures,
                thermochemistry=thermochemistry,
                equilibrium_volume=equilibrium_volume,
                helmholtz_free_energy=helmholtz_free_energy,
                gibbs_free_energy=gibbs_free_energy,
            )
    
            return result
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram.html#matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram)def draw_phase_diagram(
        labels: List[str],
        gibbs_free_energies: List[PostPhononQHAGibbsFeatureResult],
        interval: int = 501,
    ) -> plt.Figure:
        """Draw the phase diagram from the Gibbs free energy obtained from PostPhononQHAGibbsFeature.
    
        Args:
            labels (list[str]): The name of each phase.
            gibbs_free_energies (List[PostPhononQHAGibbsFeatureResult]):
                The Gibbs free energy calculation results for each phase.
            interval (int, optional):
                Number of inteploate points in both temperature and pressure.
                Defaults to 501.
        Returns:
            plt.Figure : The phase diagram.
        """
        assert len(labels) == len(gibbs_free_energies)
        n_phases = len(labels)
    
        t_min = min([min(gibbs_free_energy.temperatures) for gibbs_free_energy in gibbs_free_energies])
        t_max = max([max(gibbs_free_energy.temperatures) for gibbs_free_energy in gibbs_free_energies])
        p_min = min([min(gibbs_free_energy.pressures) for gibbs_free_energy in gibbs_free_energies])
        p_max = max([max(gibbs_free_energy.pressures) for gibbs_free_energy in gibbs_free_energies])
        t_new = np.linspace(t_min, t_max, interval, endpoint=True)
        p_new = np.linspace(p_min, p_max, interval, endpoint=True)
    
        G = np.empty([n_phases, len(t_new), len(p_new)])
        for i, gibbs_free_energy in enumerate(gibbs_free_energies):
            t = np.array(gibbs_free_energy.temperatures)
            p = np.array(gibbs_free_energy.pressures)
            g = gibbs_free_energy.gibbs_free_energy
            G[i] = _interpolate(t, p, g, t_new, p_new)
    
        ind = np.any(np.isnan(G), axis=0)
        phase_diagram = np.empty([len(t_new), len(p_new)], dtype=int)
        phase_diagram[ind] = 0
        phase_diagram[~ind] = np.argmin(G[:, ~ind], axis=0) + 1
    
        fig, ax = plt.subplots()
        cmap = matplotlib.colormaps["Set1"]
        ax.imshow(
            cmap(phase_diagram),
            origin="lower",
            vmin=0,
            vmax=1,
            cmap=cmap,
            extent=(p_min / units.GPa, p_max / units.GPa, t_min, t_max),
            aspect="auto",
        )
    
        colors = [cmap(i) for i in range(0, n_phases + 1)]
        patches = [mpatches.Patch(color=colors[i], label=l) for i, l in enumerate(["unknown"] + labels)]
    
        plt.xlabel("Pressure (GPa)")
        plt.ylabel("Temperature (K)")
        plt.legend(handles=patches)
        return fig, ax
    
    
    
    
    def _interpolate(x: np.ndarray, y: np.ndarray, z: np.ndarray, x_new: np.ndarray, y_new: np.ndarray) -> np.ndarray:
        x_, y_ = np.meshgrid(y, x)
        points = np.c_[x_[~np.isnan(z)], y_[~np.isnan(z)]]
        f = CloughTocher2DInterpolator(points, z[~np.isnan(z)])
        points_new = np.stack([i.ravel() for i in np.meshgrid(y_new, x_new)]).T
        z_new: np.ndarray = f(points_new).reshape([len(x_new), len(y_new)])
        return z_new
    
