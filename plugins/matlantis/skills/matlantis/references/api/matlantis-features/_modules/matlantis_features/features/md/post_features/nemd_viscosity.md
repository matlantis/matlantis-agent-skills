# Source code for matlantis_features.features.md.post_features.nemd_viscosity
    
    
    import warnings
    from dataclasses import dataclass
    from typing import Any, Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase import Atoms
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md.ase_integrators import LangevinIntegrator
    from matlantis_features.features.md.md import MDFeature, MDFeatureResult
    from matlantis_features.features.md.md_integrator_base import MDIntegratorBase
    from matlantis_features.features.md.md_system_base import MDSystemBase
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import viscosity_converter
    
    from .rnemd_extension import RNEMDExtension, RNEMDIntegratorError, RNEMDTrajectoryError, fit_slope
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult)@dataclass
    class PostNEMDViscosityFeatureResult:
        """A dataclass for result of PostNEMDViscosityFeature."""
    
        viscosity: Dict[str, List[float]]
        viscosity_std: Dict[str, List[float]]
        slab_velocity_x: List[float]
        slab_z: List[float]
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
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
            fig.add_trace(go.Scatter(x=self.slab_velocity_x, y=self.slab_z, showlegend=False))
            fig.update_xaxes(title="V_x sqrt(eV/amu)")
            fig.update_yaxes(title="Z Position (A)")
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Convert the result to the Pandas DataFrame object.
    
            Args:
                csv_name (str or None, optional):
                  CSV file name to write the DataFrame. If None, no CSV file is generated.
                  Defaults to None.
            Returns:
                pd.DataFrame : Resulting Pandas DataFrame object.
            """
            index = (
                ["time", "viscosity", "viscosity_std"]
                + [f"T_slab_{i}" for i in range(len(self.slab_velocity_x))]
                + [f"z_slab_{i}" for i in range(len(self.slab_z))]
            )
            data = [0.0, self.viscosity["mPa s"][0], self.viscosity_std["mPa s"][0]] + self.slab_velocity_x + self.slab_z
    
            df = pd.DataFrame({"total": data}, index=index)
            if csv_name is not None:
                df.to_csv(csv_name)
            return df
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeatureResult)@dataclass
    class ComplexNEMDViscosityFeatureResult:
        """A dataclass for result of ComplexNEMDViscosityFeature."""
    
        run_md_feature_result: MDFeatureResult
        post_nemd_viscosity_result: PostNEMDViscosityFeatureResult
    
    
    
    
    def _get_profile(atoms: Atoms) -> List[float]:
        z = atoms.get_scaled_positions()[:, 2]
        n_slab = atoms.info["n_slab"]
        slab = np.floor(z * n_slab).astype(int)
        vx = atoms.get_velocities()[:, 0]
        vx_mean = []
        for i in range(n_slab // 2 + 1):
            if (i == 0) or (i == n_slab // 2):
                vx_mean.append(vx[slab == i].mean())
            else:
                vx_mean.append((vx[slab == i].mean() + vx[slab == n_slab - i].mean()) / 2.0)
        return vx_mean
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosityFeature)class PostNEMDViscosityFeature(FeatureBase):
        """Feature class for calculating the viscosity and related properties, \
        using the non-equilibrium MD trajectory.
    
        References:
            This implementation is based on the literature:
            J. Chem. Phys. 116, 3362 (2002); https://aip.scitation.org/doi/10.1063/1.1436124
    
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/viscosity_nemd.html>`_.
        """
    
        def __call__(self, md_results: MDFeatureResult, init_time: float = 0.0) -> PostNEMDViscosityFeatureResult:
            """Calculate the viscosity from the non-equilibrium MD trajectory.
    
            Args:
                md_results (MDFeatureResult):
                  MD feature's result object containing the trajectory.
                  The simulation must be run with the rNEMD extension.
                init_time (float, optional):
                  Trajectory frames before this time is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
            Returns:
                PostNEMDViscosityFeatureResult : Result dataclass object.
            """
            traj = md_results.get_traj_obj(self.get_savedir_from_ctx())
            if "exchanged_momenta" not in traj[-1].info.keys():
                raise RNEMDTrajectoryError("No 'exchanged_momenta' property is found in the trajectory.Please run MD with the rNEMD extension")
            n_slab = traj[-1].info["n_slab"]
    
            time = []
            exchanged = []
            z = []
            vx = []
            a = []
            for atoms in traj:
                if atoms.info["time"] >= init_time:
                    assert n_slab == atoms.info["n_slab"]
                    time.append(atoms.info["time"])
                    exchanged.append(atoms.info["exchanged_momenta"])
                    vx.append(_get_profile(atoms))
                    z.append(atoms.get_cell()[2, 2])
                    a.append(atoms.get_volume() / z[-1])
    
            if len(time) == 0:
                raise ValueError("Inapproperiate value for 'init_time'.")
    
            Z = np.mean(z)  # average length alone Z-axis
            A = np.mean(a)  # average area along XY-axis
            slab_vx = np.mean(vx, axis=0)  # velocity along X-axis in each slab, unit sqrt(eV/amu)
            slab_z = (np.arange(n_slab) / n_slab + 0.5 / n_slab) * Z
            G, G_std = fit_slope(slab_z[1 : n_slab // 2 - 1], slab_vx[1 : n_slab // 2 - 1])
            J, J_std = fit_slope(2 * np.array(time) * A, np.array(exchanged))
    
            viscosity = J / G  # unit amu A^-1 fs^-1
            viscosity_std = (J_std / J + G_std / G) * viscosity
            return PostNEMDViscosityFeatureResult(
                viscosity=viscosity_converter([viscosity.tolist()]),
                viscosity_std=viscosity_converter([viscosity_std.tolist()]),
                slab_velocity_x=slab_vx.tolist(),
                slab_z=slab_z[: n_slab // 2 + 1].tolist(),
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature.html#matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature)class ComplexNEMDViscosityFeature(FeatureBase):
        """The matlantis-feature for calculating viscosity related properties.
    
        This feature calls both MDFeature with RNEMDExtension and PostNEMDViscosityFeature.
        """
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature.html#matlantis_features.features.md.post_features.nemd_viscosity.ComplexNEMDViscosityFeature.__init__)    def __init__(
            self,
            integrator: MDIntegratorBase,
            n_run: int,
            rnemd_n_slab: int = 20,
            rnemd_interval: int = 100,
            init_time: float = 0.0,
            md_feature_args: Optional[Dict[str, Any]] = None,
        ):
            """Initialize an instance.
    
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
                  Trajectory frames before this time (in fs) is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                md_feature_args (dict[str, Any] or None, optional): Additional options for MDFeature
            """
            if isinstance(integrator, LangevinIntegrator):
                raise RNEMDIntegratorError("rNEMD cannot be used with Langevin dynamics")
    
            super(ComplexNEMDViscosityFeature, self).__init__()
            if md_feature_args is None:
                mdargs = {}
            else:
                mdargs = md_feature_args
            with self.init_scope():
                self.md = MDFeature(integrator=integrator, n_run=n_run, traj_init_frame=True, **mdargs)
                self.post_viscosity = PostNEMDViscosityFeature()
                self.rnemd = RNEMDExtension(rnemd_type="viscosity", n_slab=rnemd_n_slab)
                self.rnemd_interval = rnemd_interval
                self.init_time = init_time
    
    
    
        def __call__(self, system: MDSystemBase) -> ComplexNEMDViscosityFeatureResult:
            """Run the rNEMD simulation and calculate the viscosity.
    
            Args:
                system (MDSystemBase): Target system to calculate the viscosity.
            Returns:
                ComplexNEMDViscosityFeatureResult : Result dataclass object.
            """
            md_result = self.md(system, extensions=[(self.rnemd, self.rnemd_interval)])
            viscosity_res = self.post_viscosity(md_result, init_time=self.init_time)
            return ComplexNEMDViscosityFeatureResult(run_md_feature_result=md_result, post_nemd_viscosity_result=viscosity_res)
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=0.1)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity)class PostNEMDViscosity(PostNEMDViscosityFeature):
        """The old version of PostNEMDViscosityFeature. \
        This class is deprecated. Please use PostNEMDViscosityFeature instead."""
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity.html#matlantis_features.features.md.post_features.nemd_viscosity.PostNEMDViscosity.__init__)    def __init__(self) -> None:
            super(PostNEMDViscosity, self).__init__()
            warnings.warn(
                "This class is deprecated. Please use PostEMDViscosityFeature instead.",
                FutureWarning,
            )
    
    
    
