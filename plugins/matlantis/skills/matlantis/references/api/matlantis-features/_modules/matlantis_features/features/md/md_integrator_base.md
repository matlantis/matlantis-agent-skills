# Source code for matlantis_features.features.md.md_integrator_base
    
    
    from .md_system_base import MDStateDict
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase)class MDIntegratorBase:
        """Base class for md integrator."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase.__init__)    def __init__(self, timestep: float) -> None:
            """Initialize an instance.
    
            Args:
                timestep (float): Integration time step in fs unit.
            """
            self.timestep = timestep
    
    
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the integrator.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            return dict()
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the integrator.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            pass
    
        @property
        def current_step(self) -> int:
            """Get the current integration step.
    
            Returns:
                int : Number of the current integration step.
            """
            raise NotImplementedError()
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_integrator_base.MDIntegratorBase.html#matlantis_features.features.md.md_integrator_base.MDIntegratorBase.set_temperature)    def set_temperature(self, value: float) -> None:
            """Set the target temperature of integrator.
    
            If this integrator does not support temperature control,
            this method raises the NotImplementedError exception.
    
            Args:
                value (float): Temperature in K unit.
            """
            raise NotImplementedError("this ensemble does not support temperature control")
    
    
    
