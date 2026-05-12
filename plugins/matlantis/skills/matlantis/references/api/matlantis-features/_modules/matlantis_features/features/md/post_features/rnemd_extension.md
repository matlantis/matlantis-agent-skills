# Source code for matlantis_features.features.md.post_features.rnemd_extension
    
    
    import warnings
    from typing import Tuple, cast
    
    import numpy as np
    from numpy.exceptions import RankWarning
    
    from matlantis_features.features.md import MDExtensionBase
    from matlantis_features.features.md.ase_md_system import ASEMDSystem
    from matlantis_features.features.md.md_integrator_base import MDIntegratorBase
    from matlantis_features.features.md.md_system_base import MDSystemBase
    
    
    class RNEMDError(Exception):  ## noqa
        pass
    
    
    class RNEMDTrajectoryError(RNEMDError):  ## noqa
        pass
    
    
    class RNEMDIntegratorError(RNEMDError):  ## noqa
        pass
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.html#matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension)class RNEMDExtension(MDExtensionBase):
        """Extension class for the reverse non-equilibrium MD simulation."""
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.html#matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.__init__)    def __init__(
            self,
            rnemd_type: str,
            n_slab: int,
        ):
            """Initialize an instance.
    
            Args:
                rnemd_type (str): Name of the rNEMD calculation type
                ("viscosity" or "thermal_conductivity").
                n_slab (int): Number of the slabs.
            """
            assert n_slab % 2 == 0
            assert rnemd_type in ["viscosity", "thermal_conductivity"]
            self.rnemd_type = rnemd_type
            self.n_slab = n_slab
            self.slab_size = 1.0 / self.n_slab
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.html#matlantis_features.features.md.post_features.rnemd_extension.RNEMDExtension.__call__)    def __call__(self, system: MDSystemBase, integrator: MDIntegratorBase) -> None:
            """Perform the exchange process of the rNEMD simulation.
    
            Args:
                system (MDSystemBase): Target system of the rNEMD simulation.
                integrator (MDIntegratorBase): Integrator used for the rNEMD simulation.
            """
            ase_sys = cast(ASEMDSystem, system)
            atoms = ase_sys.ase_atoms
    
            info = atoms.info
            if "n_slab" not in info:
                info["n_slab"] = self.n_slab
            if self.rnemd_type == "viscosity" and "exchanged_momenta" not in info:
                info["exchanged_momenta"] = 0.0
            if self.rnemd_type == "thermal_conductivity" and "exchanged_kinetic" not in info:
                info["exchanged_kinetic"] = 0.0
    
            z = atoms.get_scaled_positions()[:, 2]
            atn = atoms.get_atomic_numbers()
            momenta = atoms.get_momenta()
            kinetic = 0.5 * np.sum(momenta**2, axis=1) / atoms.get_masses()
    
            ind = np.where((z >= 0.0) & (z < self.slab_size))[0]
            if len(ind) == 0:
                raise ValueError("Cannot find any atoms inside slab 0. Please reduce 'n_slab' ")
    
            if self.rnemd_type == "viscosity":
                exchange_id_0 = ind[momenta[ind, 0].argmax()]
            else:
                exchange_id_0 = ind[kinetic[ind].argmax()]
    
            # Reognize exchange atomic type
            atom_type = atn[exchange_id_0]
    
            # Identify exchange atom in slab m
            ind = np.where(np.logical_and(z >= 0.5, z < 0.5 + self.slab_size, atn == atom_type))[0]
            if len(ind) == 0:
                raise ValueError(f"Cannot find any atoms inside slab {self.n_slab // 2}. Please reduce 'n_slab' ")
    
            if self.rnemd_type == "viscosity":
                exchange_id_1 = ind[momenta[ind, 0].argmin()]
            else:
                exchange_id_1 = ind[kinetic[ind].argmin()]
    
            # Accumulate exchanged kinetic energy
            if self.rnemd_type == "viscosity":
                atoms.info["exchanged_momenta"] += momenta[exchange_id_0, 0] - momenta[exchange_id_1, 0]
            else:
                atoms.info["exchanged_kinetic"] += kinetic[exchange_id_0] - kinetic[exchange_id_1]
    
            # swap
            if self.rnemd_type == "viscosity":
                momenta[[exchange_id_0, exchange_id_1], 0] = momenta[[exchange_id_1, exchange_id_0], 0]
            else:
                momenta[[exchange_id_0, exchange_id_1]] = momenta[[exchange_id_1, exchange_id_0]]
    
            atoms.set_momenta(momenta)
    
    
    
    
    
    
    [[docs]](../../../../../generated/matlantis_features.features.md.post_features.rnemd_extension.fit_slope.html#matlantis_features.features.md.post_features.rnemd_extension.fit_slope)def fit_slope(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Perform the linear fitting for the post rNEMD analysis.
    
        Args:
            x (np.ndarray): X axis data.
            y (np.ndarray): Y axis data.
        Returns:
            tuple[np.ndarray, np.ndarray] :
              Result of the linear fitting (slope and slope std).
        """
        with warnings.catch_warnings():
            warnings.simplefilter("error", category=RankWarning)
            try:
                (slope, intercept), cov = np.polyfit(x, y, 1, cov=True)
            except RankWarning:
                raise ValueError(
                    "This exception usually occurs because the simulation box is divided "
                    "into too few slabs (n_slab <= 8). "
                    "Please increase the number of slabs 'n_slab' when performing NEMD simulation. "
                )
        slope_std, intercept_std = np.sqrt(np.diag(cov))
        return slope, slope_std
    
    
    
