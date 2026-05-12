# Source code for matlantis_features.features.md.post_features.specific_heat
    
    
    import warnings
    from dataclasses import dataclass
    from typing import Dict, List, Optional
    
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.md import MDFeatureResult
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.units_util import specific_heat_converter
    
    from .thermal_expansion import _check_traj, _extract_traj
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.specific_heat.EnthalpyResult.html#matlantis_features.features.md.post_features.specific_heat.EnthalpyResult)@dataclass
    class EnthalpyResult:
        """Result dataclass for the enthalpy."""
    
        temperature: List[float]
        enthalpy: List[float]
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult.html#matlantis_features.features.md.post_features.specific_heat.SpecificHeatResult)@dataclass
    class SpecificHeatResult:
        """Result dataclass for the specific heat."""
    
        temperature: Dict[str, List[float]]
        specific_heat: Dict[str, List[float]]
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult.html#matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult)@dataclass
    class PostMDSpecificHeatFeatureResult:
        """Result dataclass for the specific heat and enthalpy."""
    
        enthalpy: EnthalpyResult
        specific_heat: SpecificHeatResult
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult.html#matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
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
    
            text = [f"MD_{i + 1}" for i in range(len(self.enthalpy.temperature))]
            fig.add_trace(
                go.Scatter(
                    x=self.enthalpy.temperature,
                    y=self.enthalpy.enthalpy,
                    text=text,
                    mode="markers+text",
                    showlegend=False,
                ),
                row=1,
                col=1,
            )
    
            fig.add_trace(
                go.Scatter(
                    x=self.specific_heat.temperature["K"],
                    y=self.specific_heat.specific_heat["J/K/g"],
                    name="specific heat",
                ),
                row=1,
                col=2,
            )
    
            fig.update_xaxes(title="temperature (K)", row=1, col=1)
            fig.update_yaxes(title="enthalpy (eV)", row=1, col=1)
    
            fig.update_xaxes(title="temperature (K)", row=1, col=2)
            fig.update_yaxes(title="specific heat (J/K/g)", row=1, col=2)
    
            fig.update_layout(height=600, width=1200)
    
            fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeature.html#matlantis_features.features.md.post_features.specific_heat.PostMDSpecificHeatFeature)class PostMDSpecificHeatFeature(FeatureBase):
        """The matlantis-feature for calculating the specific heat from \
        the MD trajectories at different temperatures."""
    
        def __call__(
            self,
            md_results_list: List[MDFeatureResult],
            temperatures: List[float],
            init_time: float = 0.0,
            fitting_order: int = 3,
        ) -> PostMDSpecificHeatFeatureResult:
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
                PostMDSpecificHeatFeatureResult : Result dataclass object.
            """
            if len(md_results_list) <= fitting_order:
                raise ValueError(
                    "Please input more 'md_results' or reduce the 'fitting_order', "
                    "The specific heat is calculated from the polynomial fitting of T-H "
                    "relationship, The 'RankingWarning' will occur if the number of input "
                    "MD trajectories are too few or the 'fitting_order' is too large."
                )
    
            savedir = self.get_savedir_from_ctx()
            _check_traj(savedir, md_results_list)
            t = _extract_traj(savedir, md_results_list, init_time, extract_fn=lambda x: x.get_temperature())
            h = _extract_traj(
                savedir,
                md_results_list,
                init_time,
                extract_fn=lambda x: (x.get_total_energy() + x.get_volume() * x.get_stress(voigt=False, include_ideal_gas=True).trace() / 3.0)
                / x.get_masses().sum(),
            )
    
            if len(t) <= fitting_order:
                # in this case EOS cannot be fitted at each temperature.
                raise ValueError("Failed get effective enthalpy from input. Please check the input md trajectries carefully.")
    
            z = np.polyfit(np.array(t), np.array(h), fitting_order)
            fit_func = np.poly1d(z)
            fit_deriv = fit_func.deriv(1)
    
            if max(temperatures) > max(t) or min(temperatures) < min(t):
                warnings.warn(
                    f"Trying to estimated the specific heat outside of the simulation temperature range {max(t):.3f} ~ {min(t):.3f}",
                    UserWarning,
                )
    
            ent_res = EnthalpyResult(
                temperature=t,
                enthalpy=h,
            )
            sh_res = SpecificHeatResult(
                temperature={"K": temperatures},
                specific_heat=specific_heat_converter(fit_deriv(temperatures).tolist()),
            )
            return PostMDSpecificHeatFeatureResult(
                enthalpy=ent_res,
                specific_heat=sh_res,
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
