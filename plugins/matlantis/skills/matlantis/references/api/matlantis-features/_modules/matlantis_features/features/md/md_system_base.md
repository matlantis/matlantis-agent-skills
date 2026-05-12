# Source code for matlantis_features.features.md.md_system_base
    
    
    from typing import Any, Dict
    
    MDStateDict = Dict[str, Any]
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_system_base.MDSystemBase.html#matlantis_features.features.md.md_system_base.MDSystemBase)class MDSystemBase:
        """Base class for the simulation system."""
    
        @property
        def current_step(self) -> int:
            """Get the current time step of this simulation.
    
            Returns:
                int : Current time step of this simulation.
            """
            raise NotImplementedError()
    
        @property
        def current_total_step(self) -> int:
            """Get the current time step of the total simulation.
    
            If the simulation is restarted from the checkpoint,
            this method returns the total time steps,
            including the simulations before the restarts.
    
            Returns:
                int : Current time step of the total simulation.
            """
            raise NotImplementedError()
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the system.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            raise NotImplementedError()
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the system.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            raise NotImplementedError()
    
    
    
