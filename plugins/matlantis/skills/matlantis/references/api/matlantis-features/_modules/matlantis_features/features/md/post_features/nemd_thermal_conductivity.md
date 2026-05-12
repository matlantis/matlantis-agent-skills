# Source code for matlantis_features.features.md.post_features.nemd_thermal_conductivity
    
    
    from dataclasses import dataclass
    from typing import Any, Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase import Atoms, units
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md.ase_integrators import LangevinIntegrator
    from matlantis_features.features.md.md import MDFeature, MDFeatureResult
    from matlantis_features.features.md.md_integrator_base import MDIntegratorBase
    from matlantis_features.features.md.md_system_base import MDSystemBase
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import thermal_conductivity_converter
    
    from .rnemd_extension import RNEMDExtension, RNEMDIntegratorError, RNEMDTrajectoryError, fit_slope
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult)@dataclass
    class PostNEMDThermalConductivityFeatureResult:
        """Result dataclass for the post-rNEMD thermal conductivity feature."""
    
        thermal_conductivity: Dict[str, List[float]]
        thermal_conductivity_std: Dict[str, List[float]]
        slab_T: List[float]
        slab_z: List[float]
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plot the result.
    
            Args:
                plt_name (str or None, optional):
                  File name to write the plot. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                  'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                  Defaults to None.
            Returns:
                go.Figure : Resulting plotly's graph object.
            """
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=self.slab_T, y=self.slab_z, showlegend=False))
            fig.update_xaxes(title="Temperature (K)")
            fig.update_yaxes(title="Z Position (A)")
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Convert the result to the Pandas DataFrame object.
    
            Args:
                csv_name (str or None, optional):
                  CSV file name to write the DataFrame. If None, no CSV file is generated.
                  Defaults to None.
            Returns:
                pd.DataFrame : Resulting Pandas DataFrame object.
            """
            index = (
                ["time", "thermal_conductivity", "thermal_conductivity_std"]
                + [f"T_slab_{i}" for i in range(len(self.slab_T))]
                + [f"z_slab_{i}" for i in range(len(self.slab_z))]
            )
            data = [0.0, self.thermal_conductivity["W/m/K"][0], self.thermal_conductivity_std["W/m/K"][0]] + self.slab_T + self.slab_z
    
            df = pd.DataFrame({"total": data}, index=index)
            if csv_name is not None:
                df.to_csv(csv_name)
            return df
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeature.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.PostNEMDThermalConductivityFeature)class PostNEMDThermalConductivityFeature(FeatureBase):
        """The matlantis-feature for calculating the thermal conductivity from \
        the non-equilibrium MD trajectory.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/thermal_conduct.html>`_.
        """
    
        def __call__(self, md_results: MDFeatureResult, init_time: float = 0.0) -> PostNEMDThermalConductivityFeatureResult:
            """Call function.
    
            Args:
                md_results (MDFeatureResult):
                  MD feature's result object containing the trajectory.
                  The simulation must be run with the rNEMD extension.
                init_time (float, optional):
                  Trajectory frames before this time is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
            Returns:
                PostNEMDThermalConductivityFeatureResult : Result dataclass object.
            """
            traj = md_results.get_traj_obj(self.get_savedir_from_ctx())
            if "exchanged_kinetic" not in traj[-1].info.keys():
                raise RNEMDTrajectoryError("No 'exchanged_kinetic' property is found in the trajectory.Please run MD with the rNEMD extension")
            n_slab = traj[-1].info["n_slab"]
    
            time = []
            exchanged = []
            z = []
            t = []
            a = []
            for atoms in traj:
                if atoms.info["time"] >= init_time:
                    assert n_slab == atoms.info["n_slab"]
                    time.append(atoms.info["time"])
                    exchanged.append(atoms.info["exchanged_kinetic"])
                    t.append(self._get_profile(atoms))
                    z.append(atoms.get_cell()[2, 2])
                    a.append(atoms.get_volume() / z[-1])
    
            if len(time) == 0:
                raise ValueError("Inapproperiate value for 'init_time'.")
    
            Z = np.mean(z)  # average length alone Z-axis
            A = np.mean(a)  # average area along XY-axis
            slab_t = np.mean(t, axis=0)  # temperature in each slab, unit K
            slab_z = (np.arange(n_slab) / n_slab + 0.5 / n_slab) * Z
            G, G_std = fit_slope(slab_z[1 : n_slab // 2 - 1], slab_t[1 : n_slab // 2 - 1])  # unit K / A
            J, J_std = fit_slope(2 * np.array(time) * A, np.array(exchanged))  # unit eV / fs / A^2
    
            thermal_conductivity = J / G  # unit eV / K / fs / A
            thermal_conductivity_std = (J_std / J + G_std / G) * thermal_conductivity
    
            results = PostNEMDThermalConductivityFeatureResult(
                thermal_conductivity=thermal_conductivity_converter([thermal_conductivity.tolist()]),
                thermal_conductivity_std=thermal_conductivity_converter([thermal_conductivity_std.tolist()]),
                slab_T=slab_t.tolist(),
                slab_z=slab_z[: n_slab // 2 + 1].tolist(),
            )
    
            return results
    
        def _get_profile(self, atoms: Atoms) -> List[float]:
            z = atoms.get_scaled_positions()[:, 2]
            n_slab = atoms.info["n_slab"]
            slab = np.floor(z * n_slab).astype(int)
            momenta = atoms.get_momenta()
            masses = atoms.get_masses()
    
            t_mean = []
            for i in range(n_slab // 2 + 1):
                if (i == 0) or (i == n_slab // 2):
                    kinetic = np.sum(momenta[slab == i] ** 2, axis=1) / masses[slab == i] / 2.0
                    t_mean.append(kinetic.mean() / 1.5 / units.kB)
                else:
                    kinetic = np.sum(momenta[slab == i] ** 2, axis=1) / masses[slab == i] / 2.0
                    kinetic_ = np.sum(momenta[slab == n_slab - i] ** 2, axis=1) / masses[slab == n_slab - i] / 2.0
                    t_mean.append((kinetic.mean() + kinetic_.mean()) / 1.5 / units.kB / 2.0)
            return t_mean
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeatureResult.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeatureResult)@dataclass
    class ComplexNEMDThermalConductivityFeatureResult:
        """Result dataclass for the complex NEMD thermal conductivity feature."""
    
        run_md_feature_result: MDFeatureResult
        post_nemd_thermal_conductivity_result: PostNEMDThermalConductivityFeatureResult
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeature.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeature)class ComplexNEMDThermalConductivityFeature(FeatureBase):
        """The matlantis-feature for calculating the thermal conductivity from \
        the non-equilibrium MD trajectory.
    
        This feature calls both MDFeature with RNEMDExteision and PostNEMDThermalConductivityFeature.
        """
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeature.html#matlantis_features.features.md.post_features.nemd_thermal_conductivity.ComplexNEMDThermalConductivityFeature.__init__)    def __init__(
            self,
            integrator: MDIntegratorBase,
            n_run: int,
            rnemd_n_slab: int = 20,
            rnemd_interval: int = 100,
            init_time: float = 0.0,
            md_feature_args: Optional[Dict[str, Any]] = None,
        ):
            """Create the feature object.
    
            Args:
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation run.
                n_run (int):
                  Number of time steps of the MD simulation.
                rnemd_n_slab (int, optional):
                  Number of the slabs used for the rNEMD calculation. Defaults to 20.
                rnemd_interval (int, optional):
                  Number of the timestep interval for the rNEMD calculation.
                  Defaults to 100.
                init_time (float, optional):
                  Trajectory frames before this time is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                md_feature_args (dict[str, Any] or None, optional): Additional options for MDFeature
            """
            if isinstance(integrator, LangevinIntegrator):
                # From Li-san's slack comment:
                #   The rNEMD method cannot works together with Langevin themostat
                #   (I noticed this issue before, but I forget it again this time.),
                #   because, I guess, the Langevin introduce the random velocity during velocity update,
                #   and such process broken the velocity profile in rNEMD.
                #   Usually, the rNEMD+Langevin gives 4~5 times larger results of viscosity.
                raise RNEMDIntegratorError("rNEMD cannot be used with Langevin dynamics")
    
            super(ComplexNEMDThermalConductivityFeature, self).__init__()
            if md_feature_args is None:
                mdargs = {}
            else:
                mdargs = md_feature_args
            with self.init_scope():
                self.md = MDFeature(integrator=integrator, n_run=n_run, traj_init_frame=True, **mdargs)
                self.post_thermal_conductivity = PostNEMDThermalConductivityFeature()
                self.rnemd = RNEMDExtension(rnemd_type="thermal_conductivity", n_slab=rnemd_n_slab)
                self.rnemd_interval = rnemd_interval
                self.init_time = init_time
    
    
    
        def __call__(self, system: MDSystemBase) -> ComplexNEMDThermalConductivityFeatureResult:
            """Run the rNEMD simulation and calculate the thermal conductivity.
    
            Args:
                system (MDSystemBase): Target system to calculate the viscosity.
            Returns:
                ComplexNEMDThermalConductivityFeatureResult : Result dataclass object.
            """
            md_result = self.md(system, extensions=[(self.rnemd, self.rnemd_interval)])
            post_result = self.post_thermal_conductivity(md_result, init_time=self.init_time)
            return ComplexNEMDThermalConductivityFeatureResult(
                run_md_feature_result=md_result, post_nemd_thermal_conductivity_result=post_result
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=0.1)
    
    
    
