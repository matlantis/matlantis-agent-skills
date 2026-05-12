# Source code for matlantis_features.features.md.ase_md_system
    
    
    import logging
    import pickle
    from typing import Optional, Union
    
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
    from numpy.random import RandomState
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    
    from .md_integrator_base import MDIntegratorBase
    from .md_system_base import MDStateDict, MDSystemBase
    
    logger = logging.getLogger(__name__)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem)class ASEMDSystem(MDSystemBase):
        """MD System using ASE backend."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem.__init__)    def __init__(self, atoms: Union[ASEAtoms, MatlantisAtoms], step: int = 0, time: float = 0.0):
            """Create MD System using ASE backend.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms):
                  ASE's Atoms object representing the simulation system.
                step (int, optional): How many MD steps have been run before the task. Defaults to 0.
                time (float, optional):
                  How long time (in fs) has it been simulated before the task. Defaults to 0.0.
            """
            if isinstance(atoms, MatlantisAtoms):
                ase_atoms = atoms.ase_atoms
            else:
                ase_atoms = atoms
            self.ase_atoms = ase_atoms.copy()
            self.ase_atoms.info["total_step"] = step
            self.ase_atoms.info["total_time"] = time
            self.ase_atoms.info["step"] = 0
            self.ase_atoms.info["time"] = 0.0
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDSystem.html#matlantis_features.features.md.ase_md_system.ASEMDSystem.init_temperature)    def init_temperature(
            self,
            temperature: float,
            stationary: bool = False,
            zero_rotation: bool = False,
            rng: Optional[RandomState] = None,
        ) -> None:
            """Initialize the velocities by the Maxwell-Boltzmann \
                distribution of the specified temperature.
    
            Args:
                temperature (float): Target temperature in K unit.
                stationary (bool, optional):
                  If True, set the center-of-mass momentum to zero.
                  Defaults to False.
                zero_rotation (bool, optional):
                  If True, set the total angular momentum to zero by counteracting rigid rotations.
                  Defaults to False.
                rng (RandomState or None, optional):
                  Random number generator to generate the distribution.
                  If None, the default generator of ASE implementation is used.
                  Defaults to None.
            """
            MaxwellBoltzmannDistribution(self.ase_atoms, temperature_K=temperature, rng=rng)
            if stationary:
                Stationary(self.ase_atoms)
            if zero_rotation:
                ZeroRotation(self.ase_atoms)
    
    
    
        @property
        def current_step(self) -> int:
            """Get the current time step of this simulation.
    
            Returns:
                int : Current time step of this simulation.
            """
            result: int = self.ase_atoms.info["step"]
            return result
    
        @property
        def current_time(self) -> float:
            """Get the current time of this simulation.
    
            Returns:
                float : Current time step of this simulation in fs unit.
            """
            result: float = self.ase_atoms.info["time"]
            return result
    
        @property
        def current_total_step(self) -> int:
            """Get the current time step of the total simulation.
    
            If the simulation is restarted from the checkpoint,
            this method returns the total time steps,
            including the simulations before the restarts.
    
            Returns:
                int : Current time step of the total simulation.
            """
            result: int = self.ase_atoms.info["total_step"]
            return result
    
        @property
        def current_total_time(self) -> float:
            """Get the current time of the total simulation.
    
            If the simulation is restarted from the checkpoint,
            this method returns the total time,
            including the simulations before the restarts.
    
            Returns:
                float : Current time step of the total simulation in fs unit.
            """
            result: float = self.ase_atoms.info["total_time"]
            return result
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the system.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            state = dict()
            state["positions"] = self.ase_atoms.get_positions()
            state["velocities"] = self.ase_atoms.get_velocities()
            state["cell"] = self.ase_atoms.get_cell()
            state["pbc"] = self.ase_atoms.get_pbc()
            state["total_time"] = self.ase_atoms.info["total_time"]
            state["total_step"] = self.ase_atoms.info["total_step"]
            return state
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the system.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            self.ase_atoms.set_positions(state["positions"])
            self.ase_atoms.set_velocities(state["velocities"])
            self.ase_atoms.set_cell(state["cell"])
            self.ase_atoms.set_pbc(state["pbc"])
            self.ase_atoms.info["total_time"] = state["total_time"]
            self.ase_atoms.info["total_step"] = state["total_step"]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDState.html#matlantis_features.features.md.ase_md_system.ASEMDState)class ASEMDState:
        """Serializable MD simulation's state class."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDState.html#matlantis_features.features.md.ase_md_system.ASEMDState.__init__)    def __init__(self, system_state: MDStateDict, integrator_state: MDStateDict):
            """Create the MD state object.
    
            Args:
                system_state (MDStateDict): System's state object.
                integrator_state (MDStateDict): Integrator's state object.
            """
            self.data = {"system": system_state, "integrator": integrator_state}
    
    
    
        @property
        def system_state(self) -> MDStateDict:
            """Get the MD system's state data.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            return self.data["system"]
    
        @property
        def integrator_state(self) -> MDStateDict:
            """Get the MD integrator's state data.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            return self.data["integrator"]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDState.html#matlantis_features.features.md.ase_md_system.ASEMDState.restore_to)    def restore_to(self, system: ASEMDSystem, integrator: MDIntegratorBase) -> None:
            """Restore the state data stored in this object \
            to the specified system and integrator.
    
            Args:
                system (ASEMDSystem): Target system object.
                integrator (MDIntegratorBase): Target integrator object.
            """
            system.state = self.data["system"]
            integrator.state = self.data["integrator"]
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_md_system.ASEMDState.html#matlantis_features.features.md.ase_md_system.ASEMDState.from_file)    @staticmethod
        def from_file(file_name: str) -> "ASEMDState":
            """Create the MD state object from file.
    
            Args:
                file_name (str): Path name of the checkpoint file to load.
            Returns:
                ASEMDState : Created the MD state object.
            """
            with open(file_name, "rb") as f:
                state: "ASEMDState" = pickle.load(f)
            return state
    
    
    
