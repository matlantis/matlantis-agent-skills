# Source code for matlantis_features.features.md.post_features.thermal_expansion
    
    
    import warnings
    from dataclasses import dataclass
    from pathlib import Path
    from typing import Any, Callable, Dict, List, Optional
    
    import numpy as np
    import plotly.graph_objects as go
    from ase import Atoms
    from plotly.subplots import make_subplots
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md import MDFeatureResult
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import thermal_expansion_converter
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.thermal_expansion.VolumeResult.html#matlantis_features.features.md.post_features.thermal_expansion.VolumeResult)@dataclass
    class VolumeResult:
        """Result dataclass for the volume result of the thermal expansion feature."""
    
        temperature: List[float]
        volume: List[float]
        length_x: List[float]
        length_y: List[float]
        length_z: List[float]
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult.html#matlantis_features.features.md.post_features.thermal_expansion.ThermalExpansionResult)@dataclass
    class ThermalExpansionResult:
        """Result dataclass for the thermal expansion result\
        of the thermal expansion feature."""
    
        temperature: Dict[str, List[float]]
        volumetric_expansion_coeffcient: Dict[str, List[float]]
        linear_expansion_coeffcient: Dict[str, List[float]]
        linear_expansion_coeffcient_x: Dict[str, List[float]]
        linear_expansion_coeffcient_y: Dict[str, List[float]]
        linear_expansion_coeffcient_z: Dict[str, List[float]]
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.html#matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult)@dataclass
    class PostMDThermalExpansionFeatureResult:
        """Result dataclass for the thermal expansion feature."""
    
        volume: VolumeResult
        thermal_expansion: ThermalExpansionResult
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.html#matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plot the result.
    
            Args:
                plt_name (str or None, optional):
                  File name to write the plot. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                  'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                  Defaults to None.
            Returns:
                go.Figure : Resulting plotly's graph object.
            """
            fig = make_subplots(rows=1, cols=2, specs=[[{}, {"secondary_y": True}]])
    
            text = [f"MD_{i + 1}" for i in range(len(self.volume.temperature))]
            fig.add_trace(
                go.Scatter(
                    x=self.volume.temperature,
                    y=self.volume.volume,
                    text=text,
                    mode="markers+text",
                    showlegend=False,
                ),
                row=1,
                col=1,
            )
    
            fig.add_trace(
                go.Scatter(
                    x=self.thermal_expansion.temperature["K"],
                    y=self.thermal_expansion.volumetric_expansion_coeffcient["1/K"],
                    name="volumetric expansion coefficient",
                ),
                row=1,
                col=2,
            )
    
            fig.add_trace(
                go.Scatter(
                    x=self.thermal_expansion.temperature["K"],
                    y=self.thermal_expansion.linear_expansion_coeffcient["1/K"],
                    name="linear expansion coefficient",
                ),
                row=1,
                col=2,
                secondary_y=True,
            )
    
            fig.update_xaxes(title="temperature (K)", row=1, col=1)
            fig.update_yaxes(title="volume (A^3)", row=1, col=1)
    
            fig.update_xaxes(title="temperature (K)", row=1, col=2)
            fig.update_yaxes(title="volumetric expansion coefficient (1/K)", row=1, col=2)
            fig.update_yaxes(title="linear expansion coefficient (1/K)", row=1, col=2, secondary_y=True)
    
            fig.update_layout(height=600, width=1200)
    
            fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    def _check_traj(
        savedir: Path,
        md_results_list: List[MDFeatureResult],
    ) -> None:
        n = len(md_results_list)
        final_structure = [md_results.get_traj_obj(savedir)[-1] for md_results in md_results_list]
        for i in range(n):
            for j in range(i + 1, n):
                if final_structure[i] == final_structure[j]:
                    raise ValueError(f"The {i}th and {j}th trajectory is identical. Each file should be different from the others.")
    
    
    def _extract_traj(
        savedir: Path,
        md_results_list: List[MDFeatureResult],
        init_time: float,
        extract_fn: Callable[[Atoms], Any],
    ) -> List[float]:
        result = []
        for md_results in md_results_list:
            traj = md_results.get_traj_obj(savedir)
            traj_res = [extract_fn(atoms) for atoms in traj if atoms.info["time"] >= init_time]
            if len(traj_res) == 0:
                raise ValueError("Length of traj is shorter than init_time. Please check the setting of init_time.")
            result.append(float(np.mean(traj_res)))
        return result
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature.html#matlantis_features.features.md.post_features.thermal_expansion.PostMDThermalExpansionFeature)class PostMDThermalExpansionFeature(FeatureBase):
        """The matlantis-feature for calculating the thermal expansion from \
        the MD trajectories at different temperatures."""
    
        def __call__(
            self,
            md_results_list: List[MDFeatureResult],
            temperatures: List[float],
            init_time: float = 0.0,
            fitting_order: int = 3,
        ) -> PostMDThermalExpansionFeatureResult:
            """Call function.
    
            Args:
                md_results_list (list[MDFeatureResult]):
                  MD feature's result object containing the trajectory.
                temperatures (list[float]):
                  List of the temperatures the simulations were performed.
                init_time (float, optional):
                  Trajectory frames before this time is ignored for the calculation.
                  Defaults to 0.0 (i.e., all frames are used for calculation).
                fitting_order (int, optional):
                  The order used for the polynomial fitting of the T-H relationship.
                  Defaults to 3.
            Returns:
                PostMDThermalExpansionFeatureResult : Result dataclass object.
            """
            if len(md_results_list) <= fitting_order:
                raise ValueError(
                    "Please input more 'md_results' or reduce the 'fitting_order', "
                    "The thermal expansion is calculated from the polynomial fitting of T-H "
                    "relationship, The 'RankingWarning' will occur if the number of input "
                    "MD trajectories are too few or the 'fitting_order' is too large."
                )
    
            savedir = self.get_savedir_from_ctx()
            _check_traj(savedir, md_results_list)
            t = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.get_temperature())
            v = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.get_volume())
            l0 = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.cell.lengths()[0])
            l1 = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.cell.lengths()[1])
            l2 = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.cell.lengths()[2])
    
            if len(t) <= fitting_order:
                # in this case EOS cannot be fitted at each temperature.
                raise ValueError("Failed get effective volume from input. Please check the input md trajectries carefully.")
    
            z = np.polyfit(np.array(t), np.array(v), fitting_order)
            z0 = np.polyfit(np.array(t), np.array(l0), fitting_order)
            z1 = np.polyfit(np.array(t), np.array(l1), fitting_order)
            z2 = np.polyfit(np.array(t), np.array(l2), fitting_order)
    
            fit_func_volume = np.poly1d(z)
            fit_func_length_0 = np.poly1d(z0)
            fit_func_length_1 = np.poly1d(z1)
            fit_func_length_2 = np.poly1d(z2)
    
            fit_deriv_volume = fit_func_volume.deriv(1)
            fit_deriv_length_0 = fit_func_length_0.deriv(1)
            fit_deriv_length_1 = fit_func_length_1.deriv(1)
            fit_deriv_length_2 = fit_func_length_2.deriv(1)
    
            if max(temperatures) > max(t) or min(temperatures) < min(t):
                warnings.warn(
                    f"Trying to estimated the thermal expansion outside of the simulation temperature range {max(t):.3f} ~ {min(t):.3f}",
                    UserWarning,
                )
    
            volumetric_expansion_coeffcient = (1 / fit_func_volume(temperatures)) * fit_deriv_volume(temperatures)
            linear_expansion_coeffcient_0 = (1 / fit_func_length_0(temperatures)) * fit_deriv_length_0(temperatures)
            linear_expansion_coeffcient_1 = (1 / fit_func_length_1(temperatures)) * fit_deriv_length_1(temperatures)
            linear_expansion_coeffcient_2 = (1 / fit_func_length_2(temperatures)) * fit_deriv_length_2(temperatures)
            linear_expansion_coeffcient = (linear_expansion_coeffcient_0 + linear_expansion_coeffcient_1 + linear_expansion_coeffcient_2) / 3.0
    
            volume_result = VolumeResult(temperature=t, volume=v, length_x=l0, length_y=l1, length_z=l2)
            thex_result = ThermalExpansionResult(
                temperature={"K": temperatures},
                volumetric_expansion_coeffcient=thermal_expansion_converter(volumetric_expansion_coeffcient.tolist()),
                linear_expansion_coeffcient=thermal_expansion_converter(linear_expansion_coeffcient.tolist()),
                linear_expansion_coeffcient_x=thermal_expansion_converter(linear_expansion_coeffcient_0.tolist()),
                linear_expansion_coeffcient_y=thermal_expansion_converter(linear_expansion_coeffcient_1.tolist()),
                linear_expansion_coeffcient_z=thermal_expansion_converter(linear_expansion_coeffcient_2.tolist()),
            )
    
            return PostMDThermalExpansionFeatureResult(volume=volume_result, thermal_expansion=thex_result)
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
