# Source code for matlantis_features.features.common.neb
    
    
    from copy import deepcopy
    from dataclasses import dataclass
    from typing import Any, Dict, List, Optional, Union
    
    import ase
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from ase.io.trajectory import Trajectory
    from packaging import version
    
    if version.parse(ase.__version__) >= version.parse("3.23.0"):
        from ase.mep import NEB
    else:
        from ase.neb import NEB
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase, FeatureBaseOutput, FeatureBaseResult
    from matlantis_features.features.common.opt import ASEOptFeature
    from matlantis_features.utils import FeatureCost
    from matlantis_features.utils.save_util import (
        get_calculator_info,
        get_user_info,
        get_version_info,
        to_serialized_dict,
    )
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureOutput.html#matlantis_features.features.common.neb.NEBFeatureOutput)@dataclass
    class NEBFeatureOutput(FeatureBaseOutput):
        """A dataclass for result output of NEBFeature."""
    
        neb_images: List[MatlantisAtoms]
        converged: bool
        steps: int
        energy_log: List[List[float]]
        force_log: List[List[float]]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureResult.html#matlantis_features.features.common.neb.NEBFeatureResult)@dataclass
    class NEBFeatureResult(FeatureBaseResult):
        """A dataclass for result of NEBFeature."""
    
        feature: "NEBFeature"
        output: "NEBFeatureOutput"
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureResult.html#matlantis_features.features.common.neb.NEBFeatureResult.plot)    def plot(self, plt_name: Optional[str] = None) -> go.Figure:
            """Plots the energy profile along the NEB pathway.
    
            Args:
                plt_name (str or None, optional): The file name to which the plot of the energy profile
                    will be saved. Supported formats include 'png' 'jpg' 'jpeg' 'webp' 'svg'
                    and 'pdf'. If 'None' is provided, the figure will not be saved. Defaults to None.
            Returns:
                go.Figure : The plot of the NEB energy profile as a plotly.graph_objects.Figure instance.
            """
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=np.arange(len(self.neb_images)), y=self.energy_log[-1], showlegend=False))
            fig.update_xaxes(title="reaction path")
            fig.update_yaxes(title="energy (eV)")
            if plt_name:
                fig.write_image(plt_name)
            return fig
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureResult.html#matlantis_features.features.common.neb.NEBFeatureResult.to_gif)    def to_gif(self, gif_name: Optional[str] = None, rotation: str = "0x,0y,0z") -> None:
            """Creates gif file to show the movement of atoms along the NEB pathway.
    
            Args:
                gif_name (str or None, optional): The file name to which the NEB pathway will be
                    saved. Defaults to None.
                rotation (str, optional): The viewing angle. Defaults to '0x,0y,0z'.
            """
            fig = plt.figure(facecolor="white")
            ax = fig.add_subplot()
            ase.io.write(
                gif_name,
                [atoms.ase_atoms for atoms in self.neb_images],
                format="gif",
                ax=ax,
                rotation=(rotation),
            )
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureResult.html#matlantis_features.features.common.neb.NEBFeatureResult.to_traj)    def to_traj(self, traj_name: str) -> None:
            """Creates ase.io.Trajectory to save the movement of atoms.
    
            Args:
                traj_name (str): The file name to which the NEB trajectory will be saved.
            """
            traj = Trajectory(traj_name, "w")
            for atoms in self.neb_images:
                traj.write(atoms.ase_atoms)
            traj.close()
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeatureResult.html#matlantis_features.features.common.neb.NEBFeatureResult.to_dataframe)    def to_dataframe(self, csv_name: Optional[str] = None) -> pd.DataFrame:
            """Saves the NEB information to pandas.DataFrame.
    
            Args:
                csv_name (str or None, optional): The file name to which the NEB information will be
                    saved. Defaults to None.
            Returns:
                pd.DataFrame :
                    The pandas.DataFrame which contains the energy and max force of each image during
                    the optimization.
            """
            n_images = len(self.output.neb_images)
            data = {"steps": list(range(self.output.steps + 1))}
            for i in range(n_images):
                data[f"energy_{i}"] = []
                data[f"max_force_{i}"] = []
    
            for e, f in zip(self.energy_log, self.force_log):
                for i in range(n_images):
                    data[f"energy_{i}"].append(e[i])  # type: ignore
                    data[f"max_force_{i}"].append(f[i])  # type: ignore
    
            df = pd.DataFrame(data)
            if csv_name:
                df.to_csv(csv_name)
    
            return df
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature)class NEBFeature(FeatureBase):
        """The matlantis-feature to run nudged elastic band (NEB) calculation.
    
        The NEB method interpolates a number of intermediate images between the given \
        initial and final states.
        Then the constrained optimization is performed by adding spring forces between \
        images to find the minimum energy path.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature.__init__)    def __init__(
            self,
            optimizer: ASEOptFeature,
            n_images: int = 5,
            k: float = 0.1,
            climb: bool = False,
            method: str = "aseneb",
            idpp: bool = False,
        ) -> None:
            """Initialize an instance.
    
            Args:
                optimizer (ASEOptFeature): The optimizer to relax each image.
                n_images (int, optional): The number of interpolated images between the given initial
                    and final NEB states. Defaults to 5.
                k (float, optional): The spring constant of the spring force. The unit is eV/Ang.
                    Defaults to 0.1.
                climb (bool, optional): Use the climbing image NEB. Defaults to False.
                method (str, optional): The scheme for the spring force. Must be in 'aseneb',
                    'improvedtangent' and 'eb'. Defaults to 'aseneb'.
                idpp (bool, optional): Interpolate the images between initial and final states with
                    IDPP method. Defaults to False.
            """
            super(NEBFeature, self).__init__()
            assert method in ["aseneb", "eb", "improvedtangent"]
            assert n_images >= 0
            assert k > 0.0
            with self.init_scope():
                self.optimizer = optimizer
                self.n_images = n_images
                self.k = k
                self.climb = climb
                self.method = method
                self.idpp = idpp
    
    
    
        def __call__(
            self,
            atoms_init: Optional[Union[ASEAtoms, MatlantisAtoms]] = None,
            atoms_end: Optional[Union[ASEAtoms, MatlantisAtoms]] = None,
            neb_images: Optional[List[Union[ASEAtoms, MatlantisAtoms]]] = None,
        ) -> NEBFeatureResult:
            """Run NEB calculation. If 'atoms_init' and 'atoms_end' are provided, the initial \
            guess of the NEB path is automatically generated by interpolation. Otherwise, the \
            initial guess of NEB path will be read from 'neb_images'.
    
            Args:
                atoms_init (ASEAtoms or MatlantisAtoms or None, optional):
                    The initial state of the NEB calculation. Defaults to None.
                atoms_end (ASEAtoms or MatlantisAtoms or None, optional):
                    The final state of the NEB calculation. Defaults to None.
                neb_images (list[ASEAtoms or MatlantisAtoms] or None, optional):
                    The initial guess of the NEB pathway. Defaults to None.
            Returns:
                NEBFeatureResult : NEB calculation.
            """
            call_parameters = {}
    
            if atoms_init and atoms_end:
                call_parameters["atoms_init"] = (
                    atoms_init.copy() if isinstance(atoms_init, MatlantisAtoms) else MatlantisAtoms(atoms_init).copy()
                )
                call_parameters["atoms_end"] = atoms_end.copy() if isinstance(atoms_end, MatlantisAtoms) else MatlantisAtoms(atoms_end).copy()
    
                ase_atoms_init = atoms_init.ase_atoms if isinstance(atoms_init, MatlantisAtoms) else atoms_init
                ase_atoms_end = atoms_end.ase_atoms if isinstance(atoms_end, MatlantisAtoms) else atoms_end
    
                assert len(ase_atoms_init) == len(ase_atoms_end)
                assert np.all(ase_atoms_init.numbers == ase_atoms_end.numbers)
                ase_neb_images = [ase_atoms_init] + [ase_atoms_init.copy() for _ in range(self.n_images)] + [ase_atoms_end]
                # don't set a Calculator, that is done inside the Optimizer
    
                neb = NEB(
                    ase_neb_images,
                    self.k,
                    self.climb,
                    method=self.method,
                    allow_shared_calculator=False,
                    parallel=True,
                )
                if self.idpp:
                    neb.interpolate(method="idpp")
                else:
                    neb.interpolate()
    
            elif neb_images:
                assert len(neb_images) == self.n_images + 2
                call_parameters["neb_images"] = [a.copy() if isinstance(a, MatlantisAtoms) else MatlantisAtoms(a).copy() for a in neb_images]  # type: ignore
    
                ase_neb_images = []
                for atoms in neb_images:
                    ase_atoms = atoms.ase_atoms if isinstance(atoms, MatlantisAtoms) else atoms
                    # don't set a Calculator, that is done inside the Optimizer
                    ase_neb_images.append(ase_atoms)
    
                neb = NEB(ase_neb_images, allow_shared_calculator=False, parallel=True)
    
            else:
                raise ValueError("Please provide either initial/end states or the neb images.")
    
            energy_log: List[List[float]] = list()
            force_log: List[List[float]] = list()
    
            def update_log() -> None:
                energy_log_: List[float] = [float(atoms.get_potential_energy()) for atoms in neb.images]
                energy_log.append(energy_log_)
                forces = [neb.images[0].get_forces()] + np.split(neb.get_forces(), neb.nimages - 2) + [neb.images[-1].get_forces()]
                force_log_: List[float] = [float(np.linalg.norm(f, axis=1).max()) for f in forces]
                force_log.append(force_log_)
    
            self.optimizer.attach(update_log, interval=1)
            results = self.optimizer(neb)
    
            output = NEBFeatureOutput(
                neb_images=[MatlantisAtoms.from_ase_atoms(atoms) for atoms in ase_neb_images],
                converged=results.converged,
                steps=results.steps,
                energy_log=energy_log,
                force_log=force_log,
            )
    
            return NEBFeatureResult(
                feature=deepcopy(self),
                call_params=call_parameters,
                output=output,
                calculator_info=get_calculator_info(self.optimizer.estimator_fn),
                version_info=get_version_info(),
                user_info=get_user_info(),
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): The number of atoms.
            Returns:
                FeatureCost :
                    A dataclass containing cost information for the feature.
            """
            # n_atoms should be
            return FeatureCost(n_pfp_run=0, time_additional=0.1)
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.neb.NEBFeature.html#matlantis_features.features.common.neb.NEBFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Serialize NEBFeature object to dict.
    
            Returns:
                dict[str, Any] : The dictionary containing the serialized NEBFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "optimizer": self.optimizer,
                "n_images": self.n_images,
                "k": self.k,
                "climb": self.climb,
                "method": self.method,
                "idpp": self.idpp,
            }
            return to_serialized_dict(res)
    
    
    
