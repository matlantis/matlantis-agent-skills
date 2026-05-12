# Source code for matlantis_features.features.md.ase_simulation
    
    
    import logging
    import pickle
    from functools import partial
    from typing import List, Optional, cast
    
    import numpy as np
    from ase import units
    from ase.io import Trajectory
    
    from matlantis_features.utils import get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    from .ase_integrators import ASEIntegrator, NPTIntegrator
    from .ase_md_system import ASEMDState, ASEMDSystem
    from .md_extensions import MDExtensionBase
    from .md_integrator_base import MDIntegratorBase
    from .md_simulation_base import MDExtensionList, MDSimulationBase
    from .md_system_base import MDSystemBase
    
    logger = logging.getLogger(__name__)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASETrajWriter.html#matlantis_features.features.md.ase_simulation.ASETrajWriter)class ASETrajWriter(MDExtensionBase):
        """ASE's trajectory writer extension."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASETrajWriter.html#matlantis_features.features.md.ase_simulation.ASETrajWriter.__init__)    def __init__(self, file_name: str, write_mode: str, props: List[str], save_init_frame: bool = False) -> None:
            """Create ASE's trajectory writer extension.
    
            Args:
                file_name (str):
                  Path name of the trajectory file.
                write_mode (str):
                  Mode string ("w" or "a", for write and append, respectively).
                props (list[str]):
                  List of property names to be written in the traj file.
                save_init_frame (bool, optional):
                  If True, the initial structure of the simulation is
                  also stored in the trajectory file. Defaults to False.
            """
            self.traj = Trajectory(
                file_name,
                write_mode,
                properties=props,
            )
            self.save_init_frame = save_init_frame
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASETrajWriter.html#matlantis_features.features.md.ase_simulation.ASETrajWriter.__call__)    def __call__(self, system: MDSystemBase, integrator: MDIntegratorBase) -> None:
            """Write the current trajectory frame.
    
            Args:
                system (MDSystemBase):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            ase_sys = cast(ASEMDSystem, system)
            atoms = ase_sys.ase_atoms
            istep = integrator.current_step
            if istep == 0 and not self.save_init_frame:
                logger.debug("skip init frame")
                return
            logger.debug(f"Traj step: {istep}({atoms.info['total_step']}) x0: {atoms.get_positions()[0]}")
            self.traj.write(atoms)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASETrajWriter.html#matlantis_features.features.md.ase_simulation.ASETrajWriter.close)    def close(self) -> None:
            """Close the underlying trajectory stream."""
            self.traj.close()
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASECheckpoint.html#matlantis_features.features.md.ase_simulation.ASECheckpoint)class ASECheckpoint(MDExtensionBase):
        """Checkpoint extension for resuming the ASE MD simulation."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASECheckpoint.html#matlantis_features.features.md.ase_simulation.ASECheckpoint.__init__)    def __init__(self, file_name: str) -> None:
            """Create checkpoint extension for resuming the ASE MD simulation.
    
            Args:
                file_name (str): Path of the checkpoint file to be created.
            """
            self.file_name = file_name
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASECheckpoint.html#matlantis_features.features.md.ase_simulation.ASECheckpoint.__call__)    def __call__(self, system: MDSystemBase, integrator: MDIntegratorBase) -> None:
            """Create the checkpoint file of the current simulation state.
    
            Args:
                system (MDSystemBase):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            istep = integrator.current_step
            if istep == 0:
                logger.debug("skip init frame")
                return
            logger.debug(f"create checkpoint step: {integrator.current_step} file: {self.file_name}")
            state = ASEMDState(system.state, integrator.state)
            with open(self.file_name, "wb") as f:
                pickle.dump(state, f)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASESimulation.html#matlantis_features.features.md.ase_simulation.ASESimulation)class ASESimulation(MDSimulationBase):
        """MD Simulation object using ASE MD backend."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASESimulation.html#matlantis_features.features.md.ase_simulation.ASESimulation.__init__)    def __init__(
            self,
            system: MDSystemBase,
            integrator: MDIntegratorBase,
            estimator_fn: Optional[EstimatorFnType] = None,
        ) -> None:
            """Create MD Simulation object using ASE MD backend.
    
            Args:
                system (MDSystemBase):
                  Target system of the MD simulation.
                  This object should be a subclass of ASEMDSystem
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
                  This object should be a subclass of ASEIntegrator.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            self.system = system
            self.integrator = integrator
            self.estimator_fn = estimator_fn
    
            atoms = cast(ASEMDSystem, system).ase_atoms
            if isinstance(self.integrator, NPTIntegrator) and not np.all(atoms.get_cell()[[1, 2, 2], [0, 0, 1]] == 0.0):
                raise ValueError(
                    "The NPTIntegrator requires the `cell` of the input structure to be an upper "
                    "triangular matrix. You can simply rotate the MatlantisAtoms to make its cell a "
                    "upper triangular matrix with the `rotate_atoms_to_upper()` method, or you can "
                    "rotate the ASEAtoms with the "
                    "`matlantis_features.utils.atoms_util.convert_atoms_to_upper()`."
                )
    
            self.init_time = atoms.info["total_time"]
            self.init_step = atoms.info["total_step"]
    
            self._extensions: MDExtensionList = []
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASESimulation.html#matlantis_features.features.md.ase_simulation.ASESimulation.append_extension)    def append_extension(self, extension: MDExtensionBase, interval: int) -> None:
            """Append extension object.
    
            Args:
                extension (MDExtensionBase):
                  Extension object appended to this simulation.
                interval (int):
                  Timestep interval to invoke this extension object.
            """
            self._extensions.append((extension, interval))
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_simulation.ASESimulation.html#matlantis_features.features.md.ase_simulation.ASESimulation.run)    def run(self, n_run_steps: int) -> None:
            """Run the MD simulation.
    
            Args:
                n_run_steps (int): Number of simulation steps to be run.
            """
            with get_calculator(estimator_fn=self.estimator_fn) as calc:
                atoms = cast(ASEMDSystem, self.system).ase_atoms
                atoms.calc = calc
                dyn = cast(ASEIntegrator, self.integrator).create_ase_dynamics(atoms)
    
                def state_updater() -> None:
                    istep = dyn.get_number_of_steps()
                    t = dyn.dt * istep / units.fs
                    logger.debug(f"Update step: {istep} time: {t}")
                    atoms.info["step"] = istep
                    atoms.info["time"] = t
                    atoms.info["total_step"] = self.init_step + istep
                    atoms.info["total_time"] = self.init_time + t
    
                dyn.attach(state_updater, interval=1)
    
                for extension, interval in self._extensions:
                    logger.debug(f"register extension: {extension}, {interval}")
                    dyn.attach(
                        partial(extension, system=self.system, integrator=self.integrator),  # type: ignore
                        interval=interval,
                    )
    
                if isinstance(self.integrator, NPTIntegrator):
                    # The ase.md.npt.NPT does not call the observers at step==0.
                    # To keep consistency with the other dynamics obj, we add the following line.
                    dyn.call_observers()
    
                dyn.run(steps=n_run_steps)
    
    
    
