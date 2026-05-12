# Source code for matlantis_features.features.md.md_simulation_base
    
    
    from typing import List, Tuple
    
    from .md_extensions import MDExtensionBase
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_simulation_base.MDSimulationBase.html#matlantis_features.features.md.md_simulation_base.MDSimulationBase)class MDSimulationBase:
        """Base class for MD simulator."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.md_simulation_base.MDSimulationBase.html#matlantis_features.features.md.md_simulation_base.MDSimulationBase.run)    def run(self, n_run_steps: int) -> None:
            """Run the MD simulation.
    
            Args:
                n_run_steps (int): Number of the simulation steps to run.
            """
            raise NotImplementedError()
    
    
    
    
    MDExtensionList = List[Tuple[MDExtensionBase, int]]
    
