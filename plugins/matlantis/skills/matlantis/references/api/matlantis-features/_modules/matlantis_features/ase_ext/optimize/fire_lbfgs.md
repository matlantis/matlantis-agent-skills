# Source code for matlantis_features.ase_ext.optimize.fire_lbfgs
    
    
    import warnings
    from math import inf
    from typing import Optional, Union
    
    import ase
    import numpy as np
    from ase import Atoms
    from packaging import version
    
    if version.parse(ase.__version__) >= version.parse("3.23.0"):
        from ase.filters import Filter
        from ase.mep import NEB
    else:
        from ase.constraints import Filter
        from ase.neb import NEB
    
    from ase.optimize import FIRE, LBFGS
    from ase.optimize.optimize import Optimizer
    
    AtomsLike = Union[Atoms, Filter, NEB]
    
    
    def _get_fmax(atoms: AtomsLike) -> float:
        """Get the maximum force.
    
        Args:
            atoms (AtomsLike): The input structure.
        Returns:
            float : The maximum force.
        """
        forces = atoms.get_forces()
        return float(np.sqrt((forces**2).sum(axis=1).max()))
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.ase_ext.optimize.FIRELBFGS.html#matlantis_features.ase_ext.optimize.FIRELBFGS)class FIRELBFGS(Optimizer):  # type: ignore
        """A optimization method that combine FIRE and LBFGS algorithm.
    
        In this method, a threshold value is define. The FIRE (LBFGS) optimization is used when
        the maximum force is larger (smaller) than the switch hold. If the force and potential energy
        increase is observed after switching to LBFGS, the optimization methid is be switch
        back to FIRE and the threshold value will be reduced by certain ratio.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.ase_ext.optimize.FIRELBFGS.html#matlantis_features.ase_ext.optimize.FIRELBFGS.__init__)    def __init__(
            self,
            atoms: AtomsLike,
            logfile: Optional[str] = "-",
            trajectory: Optional[str] = None,
            master: Optional[bool] = None,
            append_trajectory: bool = False,
            force_consistent: Optional[bool] = None,
            switch: float = 0.05,
            switch_decrease_rate: float = 0.9,
            switch_patience_fire: int = 10,
            switch_patience_lbfgs: int = 10,
            maxstep_fire: float = 0.2,
            maxstep_lbfgs: float = 0.2,
        ) -> None:
            """Initialize an instance.
    
            Args:
                atoms (AtomsLike): The input structure.
                logfile (str or None, optional):
                    Filename for log file. If '-' is used, stdout is used. Defaults to '-'.
                trajectory (str or None, optional):
                    Filename for trajectory file. If None, no trajectory file will be saved. Defaults
                    to None.
                master (bool or None, optional):
                    If True, the current rank will save files. If None, only rank 0 will save files.
                    Defaults to None.
                append_trajectory (bool, optional):
                    If True and trajectory file exists, trajectory will be
                    appended to the existing file. Defaults to False. Defaults to False.
                force_consistent (bool or None, optional):
                    Use force-consistent energy calls. Defaults to False.
                switch (float, optional):
                    The initial value of algorithm switch force threshold value. Defaults to 0.05.
                switch_decrease_rate (float, optional):
                    The ratio of threshold value decreasing. Defaults to 0.9.
                switch_patience_fire (int, optional):
                    Number of optimization steps without switching to LBFGS after maximum force lower
                    than the threshold value. Defaults to 10.
                switch_patience_lbfgs (int, optional):
                    Number of optimization steps without switching to FIRE after the energy increase is
                    observed. Defaults to 10.
                maxstep_fire (float, optional):
                    The maximum distance an atom can move per iteration when using FIRE.
                    Defaults to 0.2.
                maxstep_lbfgs (float, optional):
                    The maximum distance an atom can move per iteration when using LBFGS.
                    Defaults to 0.2.
            """
            if version.parse(ase.__version__) >= version.parse("3.23.0"):
                super().__init__(
                    atoms,
                    restart=None,  # NOTE: restart is not supported
                    logfile=logfile,
                    trajectory=trajectory,
                    master=master,
                    append_trajectory=append_trajectory,
                )
                if force_consistent is not None:
                    warnings.warn("force_consistent keyword is deprecated and will be ignored", FutureWarning)
            else:
                super().__init__(
                    atoms,
                    restart=None,  # NOTE: restart is not supported
                    logfile=logfile,
                    trajectory=trajectory,
                    master=master,
                    append_trajectory=append_trajectory,
                    force_consistent=force_consistent,
                )
            self._switch = switch
            self._switch_decrease_rate = switch_decrease_rate
            self._last_energy = inf
            self._min_energy_lbfgs = inf
            self._maxstep = maxstep_fire
            self._maxstep_lbfgs = maxstep_lbfgs
            self._switch_patience = switch_patience_fire
            self._switch_patience_lbfgs = switch_patience_lbfgs
            self._i_switch_patience = 0
            self._i_switch_patience_lbfgs = 0
            self.opt = FIRE(self.atoms, maxstep=self._maxstep, downhill_check=False)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.ase_ext.optimize.FIRELBFGS.html#matlantis_features.ase_ext.optimize.FIRELBFGS.step)    def step(self) -> None:
            """Take a single step."""
            self.opt.step()
            atoms: Atoms = self.atoms
            energy = atoms.get_potential_energy()
            fmax = _get_fmax(atoms)
            if isinstance(self.opt, FIRE):
                if fmax < self._switch:
                    self._i_switch_patience += 1
                else:
                    self._i_switch_patience = 0
                if self._i_switch_patience > self._switch_patience:
                    self._i_switch_patience = 0
                    self._switch *= self._switch_decrease_rate
                    self._min_energy_lbfgs = energy
                    self.opt = LBFGS(self.atoms, maxstep=self._maxstep_lbfgs)
            elif isinstance(self.opt, LBFGS):
                if energy > self._min_energy_lbfgs:
                    self._i_switch_patience_lbfgs += 1
                else:
                    self._min_energy_lbfgs = energy
                    self._i_switch_patience_lbfgs = 0
                if energy > self._last_energy and fmax > self._switch:
                    self._i_switch_patience_lbfgs = self._switch_patience_lbfgs + 1
                if self._i_switch_patience_lbfgs > self._switch_patience_lbfgs:
                    self._i_switch_patience_lbfgs = 0
                    self.opt = FIRE(self.atoms, maxstep=self._maxstep, downhill_check=False)
            self._last_energy = energy
    
    
    
