# Source code for matlantis_features.features.md.md_extensions
    
    
    import logging
    
    import numpy as np
    
    from .ase_integrators import NPTBerendsenIntegrator, NPTIntegrator
    from .ase_md_system import ASEMDSystem
    from .md_integrator_base import MDIntegratorBase
    
    logger = logging.getLogger(__name__)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase)class MDExtensionBase:
        """Base class for the MD extension."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.MDExtensionBase.html#matlantis_features.features.md.md_extensions.MDExtensionBase.__call__)    def __call__(self, system: ASEMDSystem, integrator: MDIntegratorBase) -> None:
            """Call the MD extension.
    
            Args:
                system (ASEMDSystem):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            raise NotImplementedError()
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.TemperatureScheduler.html#matlantis_features.features.md.md_extensions.TemperatureScheduler)class TemperatureScheduler(MDExtensionBase):
        """Linear temperature scheduler exteisnon.
    
        This extension linearly change the temperature of the system
        from the start to the end temperature,
        during the given time steps.
        After the given time step, the temperature is kept to the end temperature.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.TemperatureScheduler.html#matlantis_features.features.md.md_extensions.TemperatureScheduler.__init__)    def __init__(self, start_value: float, end_value: float, num_total_steps: int):
            """Create linear temperature scheduler exteisnon.
    
            Args:
                start_value (float):
                  Start temperature in K unit.
                end_value (float):
                  End temperature in K unit.
                num_total_steps (int):
                  Number of steps for the linear temperature change.
            """
            self.start_value = start_value
            self.end_value = end_value
            self.nsteps = num_total_steps
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.TemperatureScheduler.html#matlantis_features.features.md.md_extensions.TemperatureScheduler.__call__)    def __call__(self, system: ASEMDSystem, integrator: MDIntegratorBase) -> None:
            """Change the temperature of the system.
    
            Args:
                system (ASEMDSystem):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            istep = system.current_total_step
            ratio = max(0.0, min(1.0, istep / self.nsteps))
            cur_temp = (self.end_value - self.start_value) * ratio + self.start_value
            logger.debug(f"change temperature to {cur_temp}")
            integrator.set_temperature(cur_temp)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.DeformScheduler.html#matlantis_features.features.md.md_extensions.DeformScheduler)class DeformScheduler(MDExtensionBase):
        """Linear cell shape scheduler extension."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.DeformScheduler.html#matlantis_features.features.md.md_extensions.DeformScheduler.__init__)    def __init__(self, final: np.ndarray, num_total_steps: int):
            """Create linear cell shape scheduler.
    
            Args:
                final (np.ndarray): The final cell shape after deformation.
                num_total_steps (int): Number of MD steps for the cell shape change procedure.
            """
            self.final = final
            self.num_total_steps = num_total_steps
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_extensions.DeformScheduler.html#matlantis_features.features.md.md_extensions.DeformScheduler.__call__)    def __call__(self, system: ASEMDSystem, integrator: MDIntegratorBase) -> None:
            """Change the cell of the system.
    
            Args:
                system (ASEMDSystem):
                  Target system of the MD simulation.
                integrator (MDIntegratorBase):
                  Integrator used for the MD simulation.
            """
            if not hasattr(self, "start"):
                self.start = system.ase_atoms.get_cell()[:]
                self.start_step = system.current_total_step
                self.final_step = system.current_total_step + self.num_total_steps
                if isinstance(integrator, NPTIntegrator) or isinstance(integrator, NPTBerendsenIntegrator):
                    raise ValueError("DeformSchedular cannot be used in NPT MD.")
            istep = system.current_total_step
            if istep <= self.final_step:
                cell = self.start + (self.final - self.start) * (istep - self.start_step) / self.num_total_steps
                system.ase_atoms.set_cell(cell, scale_atoms=True)
    
    
    
