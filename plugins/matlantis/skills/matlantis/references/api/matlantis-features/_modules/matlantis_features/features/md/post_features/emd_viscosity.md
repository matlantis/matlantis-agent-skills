# Source code for matlantis_features.features.md.post_features.emd_viscosity
    
    
    import warnings
    from dataclasses import dataclass
    from typing import Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase import units
    from plotly.subplots import make_subplots
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md.md import MDFeatureResult
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import viscosity_converter
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult)@dataclass
    class PostEMDViscosityFeatureResult:
        """A dataclass for result of PostEMDViscosityFeature."""
    
        num_of_segments: int
        viscosity: Dict[str, List[float]]
        viscosity_std: Dict[str, List[float]]
        time: np.ndarray
        acf_segment: np.ndarray
        acf_avg: np.ndarray
        vis_segment: np.ndarray
        vis_avg: np.ndarray
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plot the result.
    
            Args:
                plt_name (str or None, optional):
                  File name to write the plot. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                  'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                  Defaults to None.
            Returns:
                go.Figure : Resulting plotly's graph object.
            """
            fig = make_subplots(rows=1, cols=2)
            for i in range(self.acf_segment.shape[1]):
                fig.add_trace(
                    go.Scatter(
                        x=self.time,
                        y=self.acf_segment[:, i],
                        mode="lines",
                        showlegend=False,
                        line=dict(color="#7f7f7f"),
                    ),
                    row=1,
                    col=1,
                )
            fig.add_trace(
                go.Scatter(x=self.time, y=self.acf_avg, showlegend=False),
                row=1,
                col=1,
            )
            for i in range(self.acf_segment.shape[1]):
                fig.add_trace(
                    go.Scatter(
                        x=self.time,
                        y=self.vis_segment[:, i] * 1e28 / units.kg,
                        mode="lines",
                        showlegend=False,
                        line=dict(color="#7f7f7f"),
                    ),
                    row=1,
                    col=2,
                )
            fig.add_trace(
                go.Scatter(x=self.time, y=self.vis_avg * 1e28 / units.kg, showlegend=False),
                row=1,
                col=2,
            )
            idx = len(self.vis_segment) // 2
            fig.add_trace(
                go.Scatter(
                    x=[self.time[idx]],
                    y=self.viscosity["mPa s"],
                    showlegend=False,
                    error_y=dict(type="data", array=self.viscosity_std["mPa s"], visible=True),
                ),
                row=1,
                col=2,
            )
            fig.update_xaxes(title="time [fs]", row=1, col=1)
            fig.update_yaxes(title="autocorrelation function", row=1, col=1)
            fig.update_xaxes(title="time [fs]", row=1, col=2)
            fig.update_yaxes(title="viscosity [mPa s]", row=1, col=2)
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Convert the result to the Pandas DataFrame object.
    
            Args:
                csv_name (str or None, optional):
                  CSV file name to write the DataFrame. If None, no CSV file is generated.
                  Defaults to None.
            Returns:
                pd.DataFrame : Resulting Pandas DataFrame object.
            """
            df = pd.DataFrame()
            df["time"] = self.time
            for i in range(self.num_of_segments):
                df[f"acf_segment_{i}"] = self.acf_segment[:, i]
            df["acf"] = self.acf_avg
            for i in range(self.num_of_segments):
                df[f"viscosity_segment_{i}"] = self.vis_segment[:, i]
            df["viscosity"] = self.vis_avg
            # current implement ommited the function fitting of viscosity
            df["viscosity_fitted"] = self.vis_avg
            if csv_name is not None:
                df.to_csv(csv_name)
            return df
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosityFeature)class PostEMDViscosityFeature(FeatureBase):
        """The matlantis-feature for calculating the viscosity and related properties, \
        using the equilibrium MD trajectory.
    
        References:
            This implementationis based on the literature:
            J. Chem. Phys. 106, 9327 (1997); https://aip.scitation.org/doi/10.1063/1.474002
    
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/viscosity_emd.html>`_.
        """
    
        def __call__(
            self,
            md_results: MDFeatureResult,
            init_time: float = 0.0,
            number_of_segments: int = 1,
        ) -> PostEMDViscosityFeatureResult:
            """Calculate the viscosity from the MD trajectory.
    
            Args:
                md_results (MDFeatureResult):
                  MD feature's result object containing the trajectory.
                init_time (float, optional):
                  Trajectory frames before this time is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                number_of_segments (int, optional):
                  Divides the given trajectory in to segments to allow statistical analysis.
                  Defaults to 1.
            Returns:
                PostEMDViscosityFeatureResult : Result dataclass object.
            """
            traj = md_results.get_traj_obj(self.get_savedir_from_ctx())
            frame_step = traj[1].info["time"] - traj[0].info["time"]
            if traj[-1].info["time"] < init_time:
                raise ValueError("init_time is larger than the length of trajectory.")
    
            # Get sliced trajectory
            for i, atoms in enumerate(traj):
                if atoms.info["time"] >= init_time:
                    st = i
                    break
    
            sliced_traj = traj[st:]
            n_images = len(sliced_traj)
            segment = n_images // number_of_segments
    
            time = np.empty(n_images)
            P = np.empty([n_images, 6])
            V = np.empty(n_images)
            T = np.empty(n_images)
            init_time = sliced_traj[0].info["time"]
            for i, atoms in enumerate(sliced_traj):
                time[i] = atoms.info["time"] - init_time
                stress = atoms.get_stress(voigt=False, include_ideal_gas=True)
                stress_ = (stress + stress.T) / 2.0 - stress.trace() / 3.0 * np.eye(3)
                P[i] = stress_[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]
                V[i] = atoms.get_volume()
                T[i] = atoms.get_temperature()
    
            # calculate viscosity in each segment
            acf_segment = np.empty([segment - 1, number_of_segments])
            vis_segment = np.empty([segment - 1, number_of_segments])
            for i in range(number_of_segments):
                acf_12 = _get_acf(P[i * segment : (i + 1) * segment, 3])
                acf_13 = _get_acf(P[i * segment : (i + 1) * segment, 4])
                acf_23 = _get_acf(P[i * segment : (i + 1) * segment, 5])
                acf_segment[:, i] = (acf_12 + acf_13 + acf_23) / 3.0  # eV^2 / A^6
                scale = V[i * segment : (i + 1) * segment].mean() / units.kB / T[i * segment : (i + 1) * segment].mean()
                vis_segment[:, i] = _trap_cum(acf_segment[:, i] * scale, frame_step) * (
                    units.kg / units.J * 1e-10
                )  # unit eV fs A^-3 -> amu/A/fs
    
            idx = len(vis_segment) // 2
            return PostEMDViscosityFeatureResult(
                num_of_segments=number_of_segments,
                viscosity=viscosity_converter([float(vis_segment[idx].mean())]),
                viscosity_std=viscosity_converter([float(vis_segment[idx].std())]),
                time=time[: segment - 1],
                acf_segment=acf_segment,
                acf_avg=acf_segment.mean(axis=1),
                vis_segment=vis_segment,
                vis_avg=vis_segment.mean(axis=1),
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
    
    def _get_acf(x: np.ndarray) -> np.ndarray:
        # Calculate the auto correlation function.
        N = len(x)
        acf = np.empty(N - 1)
    
        acf[0] = np.sum(x**2) / N
        for i in range(1, N - 1):
            acf[i] = np.sum(x[:-i] * x[i:]) / (N - i)
    
        return acf  # type: ignore
    
    
    def _trap_cum(x: np.ndarray, dx: float) -> np.ndarray:
        # Calculate the cumulative of the integral using trapezoidal rule.
        N = len(x)
        cum = np.empty(N)
    
        for i in range(N):
            cum[i] = np.trapz(x[:i], dx=dx)
    
        return cum  # type: ignore
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity)class PostEMDViscosity(PostEMDViscosityFeature):
        """The old version of PostEMDViscosityFeature. \
        This class is deprecated. Please use PostEMDViscosityFeature instead."""
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity.html#matlantis_features.features.md.post_features.emd_viscosity.PostEMDViscosity.__init__)    def __init__(self) -> None:
            super(PostEMDViscosity, self).__init__()
            warnings.warn(
                "This class is deprecated. Please use PostEMDViscosityFeature instead.",
                FutureWarning,
            )
    
    
    
