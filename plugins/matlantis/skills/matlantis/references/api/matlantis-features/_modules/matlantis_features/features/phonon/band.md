# Source code for matlantis_features.features.phonon.band
    
    
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly
    import plotly.graph_objects as go
    import plotly.subplots
    
    from matlantis_features.atoms import MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.dos import PostPhononDOSFeatureResult
    from matlantis_features.features.phonon.force_constant import ForceConstantFeatureResult
    from matlantis_features.features.phonon.utils import KptPath, PhononFrequency
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.visual_utils.brillouin_zone import get_brillouin_zone_3d
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult)@dataclass
    class PostPhononBandFeatureResult:
        """A dataclass for result of PostPhononBandFeature."""
    
        kpts: np.ndarray
        coords: np.ndarray
        special_kpt_indices: np.ndarray
        labels: List[str]
        frequency: np.ndarray
        unit: str
        unit_cell_atoms: MatlantisAtoms
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plots the phonon band structure.
    
            Args:
                plt_name (str or None, optional): The file name to which the plot of the phonon
                    band structure will be saved. Supported formats include 'png' 'jpg' 'jpeg'
                    'webp' 'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                    Defaults to None.
            Returns:
                go.Figure :
                    The plot of the phonon band structure as a plotly.graph_objects.Figure instance.
            """
            unit = self.unit
            fig = go.Figure()
            for values in self.frequency.T:
                fig.add_trace(
                    go.Scatter(
                        x=self.coords,
                        y=values,
                        mode="lines",
                        marker=dict(color=plotly.colors.qualitative.Plotly[0]),
                        showlegend=False,
                    )
                )
    
            special_coords = self.coords[self.special_kpt_indices]
            tickvals = [special_coords[0]]
            ticktext = [self.labels[0]]
            vline = []
            for i in range(1, len(self.labels)):
                if special_coords[i] == special_coords[i - 1]:
                    ticktext[-1] = f"{self.labels[i - 1]}|{self.labels[i]}"
                    vline.append(special_coords[i])
                else:
                    tickvals.append(special_coords[i])
                    ticktext.append(self.labels[i])
    
            for x in vline:
                fig.add_vline(x=x, line_width=1, line_dash="dash", line_color="gray")
    
            fig.update_layout(xaxis=dict(tickvals=tickvals, ticktext=ticktext))
            fig.update_yaxes(title=f"Frequency ({unit})")
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Converts the phonon band structure to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The csv file path to which the dataframe will be
                    saved. If 'None' is provided, the dataframe will not be saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains the coordinates of the band path, the labels
                    of special points and the phonon vibration frequencies.
            """
            n_kpts = len(self.kpts)
            labels: List[Any] = [None for _ in range(n_kpts)]
            for i, l in zip(self.special_kpt_indices, self.labels):
                labels[i] = l
    
            df = pd.DataFrame(
                {
                    "kx": self.kpts[:, 0],
                    "ky": self.kpts[:, 1],
                    "kz": self.kpts[:, 2],
                    "coords": self.coords,
                    "labels": labels,
                }
            )
            for i, values in enumerate(self.frequency.T):
                df[str(i)] = values
    
            if csv_name:
                df.to_csv(csv_name)
    
            return df
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeatureResult.html#matlantis_features.features.phonon.band.PostPhononBandFeatureResult.visual_k_path)    def visual_k_path(
            self,
            show_reciprocal_basis: bool = False,
            show_special_k_points: bool = True,
            marker_size: int = 2,
        ) -> go.Figure:
            """To visualize the k point path in the 3D plot of the first brillouin zone.
    
            Args:
                show_reciprocal_basis (bool, optional):
                  Show the basis vector of the reciprocal cell or not. Defaults to False.
                show_special_k_points (bool, optional):
                  Show the special k points or not. Defaults to True.
                marker_size (int, optional):
                  The marker size of the k points. Defaults to 2.
            Returns:
                go.Figure : The 3D plot of the k point path.
            """
            ase_atoms = self.unit_cell_atoms.ase_atoms
            brillouin_zone = get_brillouin_zone_3d(ase_atoms.cell, show_reciprocal_basis=show_reciprocal_basis)
            kpts = self.kpts.dot(ase_atoms.cell.reciprocal()[:])
            brillouin_zone.add_trace(
                go.Scatter3d(
                    x=kpts[:, 0],
                    y=kpts[:, 1],
                    z=kpts[:, 2],
                    mode="markers",
                    marker=dict(size=marker_size, color="blue"),
                    showlegend=False,
                )
            )
            if show_special_k_points:
                special_kpts = self.kpts[self.special_kpt_indices].dot(ase_atoms.cell.reciprocal()[:])
                brillouin_zone.add_trace(
                    go.Scatter3d(
                        x=special_kpts[:, 0],
                        y=special_kpts[:, 1],
                        z=special_kpts[:, 2],
                        mode="markers",
                        marker=dict(size=5, color="green"),
                        text=self.labels,
                        showlegend=False,
                    )
                )
            return brillouin_zone
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the PostPhononBandFeatureResult.
    
            Returns:
                dict[str, Any] : Phonon band structure in dictionary format.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeature.html#matlantis_features.features.phonon.band.PostPhononBandFeature)class PostPhononBandFeature(FeatureBase):
        """The matlantis-feature for calculating the phonon band structure.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/phonon_band.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.PostPhononBandFeature.html#matlantis_features.features.phonon.band.PostPhononBandFeature.__init__)    def __init__(self) -> None:
            """Initiate an instance."""
            super(PostPhononBandFeature, self).__init__()
    
    
    
        def __call__(
            self,
            force_constant: ForceConstantFeatureResult,
            labels: Optional[List[str]] = None,
            special_kpts: Optional[np.ndarray] = None,
            n_kpts: Optional[List[int]] = None,
            total_n_kpts: int = 100,
            unit: str = "meV",
            style: str = "ase",
        ) -> PostPhononBandFeatureResult:
            """Calculates the phonon band structure from the force constant.
    
            Args:
                force_constant (ForceConstantFeatureResult): The calculation result of
                    the ForceConstantFeature.
                labels (list[str] or None, optional): The labels of special points that form the band
                    path. Please use ',' '|' or '/' to represent a discontinuous jump. If None is
                    provided, the band path will be automatically generated according to the type of
                    bravais lattice. Defaults to None.
                special_kpts (np.ndarray or None, optional): The coordinates of special points in
                    reciprocal space. The 'special_kpts' should be a numpy.ndarray in the shape of
                    (N, 3), where N is the same as the length of 'labels'. If 'None' is provided,
                    the automatically generated special points will be used. Defaults to None.
                n_kpts (list[int] or None, optional): The number of interpolated points in each
                    segment of the path. If None is provided, the values will be estimated from
                    'total_n_kpts'. Defaults to None.
                total_n_kpts (int, optional): The number of interpolated points along the whole band
                    path. This parameter will be ignored if 'n_kpts' is specified. Defaults to 100.
                unit (str, optional): The unit of phonon frequency. Must be one of: 'meV', 'eV', 'THz'
                    or 'cm^-1'. Defaults to 'meV'.
                style (str, optional): If 'labels' and 'special_kpts' are all not provided,
                    the default k-point path will be used. The parameter 'style' can control which
                    kind of default path to use. Currently, 'ase', which is ASE style, and 'phonon_db',
                    which follows the Phonon database (http://phonondb.mtl.kyoto-u.ac.jp/), are
                    supported. Defaults to "ase".
            Returns:
                PostPhononBandFeatureResult : The calculation results of the phonon band structure.
            """
            pf = PhononFrequency(
                force_constant=force_constant.force_constant,
                supercell=force_constant.supercell,
                unit_cell_atoms=force_constant.unit_cell_atoms,
            )
    
            kpt_path = KptPath(pf.cell, labels, special_kpts, n_kpts, total_n_kpts, style)
            frequency = pf.get_frequency(kpt_path.kpts, unit=unit)
    
            result = PostPhononBandFeatureResult(
                **{  # type: ignore
                    "kpts": kpt_path.kpts,
                    "coords": kpt_path.coords,
                    "special_kpt_indices": kpt_path.special_kpt_indices,
                    "labels": kpt_path.labels,
                    "frequency": frequency,
                    "unit": unit,
                    "unit_cell_atoms": force_constant.unit_cell_atoms,
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
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.band.plot_band_dos.html#matlantis_features.features.phonon.band.plot_band_dos)def plot_band_dos(
        band_results: PostPhononBandFeatureResult,
        dos_results: PostPhononDOSFeatureResult,
        plt_name: Optional[str] = None,
    ) -> go.Figure:
        """Plot the phonon band structure and phonon density of states in one figure.
    
        Args:
            band_results (PostPhononBandFeatureResult): The result of PostPhononBandFeature.
            dos_results (PostPhononDOSFeatureResult): The result of PostPhononDOSFeature.
            plt_name (str or None, optional): The file name to which the figure will be saved.
                    Supported formats include 'png' 'jpg' 'jpeg' 'webp' 'svg' and 'pdf'. If
                    'None' is provided, the figure will not be saved. Defaults to None.
        Returns:
            go.Figure :
                The plot of the phonon band structure and phonon density of states as a
                plotly.graph_objects.Figure instance.
        """
        assert band_results.unit == dos_results.unit
        energy = dos_results.energy[dos_results.unit]
        y_range = [min(energy), max(energy)]
    
        band_plot = band_results.plot()
        dos_plot = dos_results.plot(reverse_axis=True)
    
        fig = plotly.subplots.make_subplots(1, 2, column_widths=[0.7, 0.3])
        for trace in band_plot.select_traces():
            fig.add_trace(trace, row=1, col=1)
    
        for shape in band_plot.select_shapes():
            fig.add_shape(shape)
    
        for trace in dos_plot.select_traces():
            fig.add_trace(trace, row=1, col=2)
    
        fig.update_xaxes(band_plot.layout.xaxis, row=1, col=1)
        fig.update_yaxes(band_plot.layout.yaxis, range=y_range, row=1, col=1)
        fig.update_xaxes(dos_plot.layout.xaxis, row=1, col=2)
        fig.update_yaxes(dos_plot.layout.yaxis, range=y_range, title="", row=1, col=2)
    
        if plt_name:
            fig.write_image(plt_name)
    
        return fig
    
    
    
