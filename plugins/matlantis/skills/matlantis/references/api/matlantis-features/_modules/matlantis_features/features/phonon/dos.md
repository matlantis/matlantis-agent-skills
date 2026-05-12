# Source code for matlantis_features.features.phonon.dos
    
    
    from dataclasses import asdict, dataclass
    from typing import Any, Dict, List, Optional
    
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase.data import chemical_symbols
    
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.phonon.force_constant import ForceConstantFeatureResult
    from matlantis_features.features.phonon.utils import PhononFrequency, get_k_mesh
    from matlantis_features.utils import FeatureCost
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult)@dataclass
    class PostPhononDOSFeatureResult:
        """A dataclass for result of PostPhononDOSFeature."""
    
        energy: Dict[str, List[float]]
        dos: Dict[str, List[float]]
        unit: str
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None, reverse_axis: bool = False) -> go.Figure:
            """Plots the phonon density of states (DOS).
    
            Args:
                plt_name (str or None, optional): The file name to which the plot of phonon
                    DOS will be saved. Supported formats include 'png' 'jpg' 'jpeg' 'webp'
                    'svg' and 'pdf'. If 'None' is provided, the figure will not be saved.
                    Defaults to None.
                reverse_axis (bool, optional): If 'False', the x-axis will be 'energy' and the y-axis
                    will be density of states. If 'True', the x-axis and y-axis will be exchanged.
    
                    Defaults to False.
    
            Returns:
                go.Figure : The plot of phonon DOS as a plotly.graph_objects.Figure instance.
            """
            keys = self.dos.keys()
            partial = len(keys) != 1
    
            fig = go.Figure()
            if partial:
                fig.add_trace(
                    go.Scatter(
                        x=self.energy[self.unit],
                        y=self.dos["total"],
                        mode="lines+markers",
                        name="total",
                    )
                )
                for key in keys:
                    if key != "total":
                        fig.add_trace(
                            go.Scatter(
                                x=self.energy[self.unit],
                                y=self.dos[key],
                                mode="lines+markers",
                                name=key,
                            )
                        )
            else:
                fig.add_trace(
                    go.Scatter(
                        x=self.energy[self.unit],
                        y=self.dos["total"],
                        mode="lines+markers",
                        showlegend=False,
                    )
                )
    
            fig.update_yaxes(title="DOS", showticklabels=False)
            fig.update_xaxes(title=f"energy ({self.unit})")
    
            if reverse_axis:
                for trace in fig.select_traces():
                    x = trace.x
                    trace.x = trace.y
                    trace.y = x
                fig.update_xaxes(title="DOS", showticklabels=False)
                fig.update_yaxes(title=f"energy ({self.unit})")
    
            if plt_name:
                fig.write_image(plt_name)
    
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.html#matlantis_features.features.phonon.dos.PostPhononDOSFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Convert the phonon DOS to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The csv file path to which the dataframe will be
                    saved. If 'None' is provided, the dataframe will not be saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains the phonon vibration energy and corresponding
                    density of states. The element projected DOS will also be saved if the system
                    contains more than one element.
            """
            df = pd.DataFrame({"energy": self.energy[self.unit]})
            for key in self.dos.keys():
                df[key] = self.dos[key]
    
            if csv_name:
                df.to_csv(csv_name, index=False)
    
            return df
    
    
    
        @property
        def dict(self) -> Dict[str, Any]:
            """Dictionary representation of the PostPhononDOSFeatureResult.
    
            Returns:
                dict[str, Any] : Phonon DOS in dictionary format.
            """
            return asdict(self)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.dos.PostPhononDOSFeature.html#matlantis_features.features.phonon.dos.PostPhononDOSFeature)class PostPhononDOSFeature(FeatureBase):
        """The matlantis-feature for calculating the phonon DOS from the force constant.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/phonon_dos.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.dos.PostPhononDOSFeature.html#matlantis_features.features.phonon.dos.PostPhononDOSFeature.__init__)    def __init__(self, partial: bool = False) -> None:
            """Initialize an instance.
    
            Args:
                partial (bool, optional): Whether to calculate the element projected DOS. Defaults to
                    False.
            """
            super(PostPhononDOSFeature, self).__init__()
            with self.init_scope():
                self.partial = partial
    
    
    
        def __call__(
            self,
            force_constant: ForceConstantFeatureResult,
            kpts: List[int],
            scheme: str = "mp",
            shift: Optional[np.ndarray] = None,
            freq_min: Optional[float] = None,
            freq_max: Optional[float] = None,
            freq_bin: Optional[float] = None,
            unit: str = "meV",
        ) -> PostPhononDOSFeatureResult:
            """Calculates the phonon density of states (DOS) from the force constant.
    
            Args:
                force_constant (ForceConstantFeatureResult): The calculation result of
                    the ForceConstantFeature.
                kpts (list[int]): The number of mesh points along each axis.
                scheme (str, optional): The type of mesh grid: either 'mp' or 'gamma'. The mesh grid is
                    determined by the Monkhorst-Pack scheme ('mp') or Gamma-centered scheme ('gamma').
                    Defaults to 'mp'.
                shift (np.ndarray or None, optional): The shift of the mesh grid in the direction along
                    the corresponding reciprocal axes. Defaults to None.
                freq_min (float or None, optional): The minimum frequency of phonon DOS calculation
                    range. Defaults to None.
                freq_max (float or None, optional): The maximum frequency of phonon DOS calculation
                    range. Defaults to None.
                freq_bin (float or None, optional): The small interval in which the density of states
                    will be calculated. Defaults to None.
                unit (str, optional): The unit of phonon frequency. Must be one of: 'meV', 'eV', 'THz'
                    or 'cm^-1'. Defaults to 'meV'.
            Returns:
                PostPhononDOSFeatureResult : The calculation results of phonon DOS.
            """
            assert scheme in ["mp", "gamma"]
            assert unit in ["meV", "eV", "THz", "cm^-1"]
            if shift is not None:
                assert shift.shape == (3,)
    
            pf = PhononFrequency(
                force_constant=force_constant.force_constant,
                supercell=force_constant.supercell,
                unit_cell_atoms=force_constant.unit_cell_atoms,
            )
            n_atoms = pf.n_atoms
            atn = pf.numbers
            spec = np.unique(atn)
            if len(spec) == 1:
                # if only has one kind of element, partial DOS is same as DOS
                partial = False
            else:
                partial = self.partial
    
            kpts = get_k_mesh(kpts, scheme, shift)  # type: ignore
            if partial:
                frequency, mode = pf.get_frequency_mode(kpts, unit=unit)  # type: ignore
            else:
                frequency = pf.get_frequency(kpts, unit=unit)  # type: ignore
    
            # range
            if freq_bin is None:
                if unit == "meV":
                    freq_bin = 1.0
                elif unit == "eV":
                    freq_bin = 0.001
                elif unit == "THz":
                    freq_bin = 0.2
                else:
                    freq_bin = 10.0
            if freq_max is None:
                freq_max = (np.floor(np.max(frequency) / freq_bin) + 2) * freq_bin
            if freq_min is None:
                freq_min = (np.floor(np.min(frequency) / freq_bin) - 1) * freq_bin
    
            if freq_min >= freq_max:
                raise ValueError(f"Please check the setting of 'freq_min' ({freq_min}) and 'freq_max' ({freq_max}). ")
    
            dos_, bins = np.histogram(
                frequency,
                range=(freq_min, freq_max),
                bins=int(np.rint((freq_max - freq_min) // freq_bin)),
            )
            energy = (bins[1:] + bins[:-1]) / 2.0
            dos = {"total": dos_ / 3 / n_atoms / len(kpts)}
    
            if partial:
                # cal weight
                weight = np.empty([frequency.shape[0], frequency.shape[1], len(spec)])
                for i in range(frequency.shape[0]):
                    for j in range(frequency.shape[1]):
                        for k, s in enumerate(spec):
                            weight[i, j, k] = np.sum(mode[i, j][atn == s] ** 2)
    
                weight /= weight.sum(axis=2, keepdims=True)
                for k, s in enumerate(spec):
                    dos_, bins = np.histogram(
                        frequency,
                        range=(freq_min, freq_max),
                        bins=int(np.rint((freq_max - freq_min) // freq_bin)),
                        weights=weight[:, :, k],
                    )
                    dos[chemical_symbols[s]] = dos_ / 3 / n_atoms / len(kpts)
    
            result = PostPhononDOSFeatureResult(energy={unit: energy.tolist()}, dos=dos, unit=unit)
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
    
    
    
