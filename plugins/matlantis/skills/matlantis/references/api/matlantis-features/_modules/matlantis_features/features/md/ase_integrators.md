# Source code for matlantis_features.features.md.ase_integrators
    
    
    import logging
    from typing import Any, Callable, Dict, Optional
    
    import numpy as np
    from ase import Atoms, units
    from ase.md.andersen import Andersen
    from ase.md.langevin import Langevin
    from ase.md.npt import NPT
    from ase.md.nptberendsen import NPTBerendsen
    from ase.md.nvtberendsen import NVTBerendsen
    from ase.md.verlet import VelocityVerlet
    from ase.optimize.optimize import Dynamics
    
    from matlantis_features.filters import IsotropicFilter
    
    from .md_integrator_base import MDIntegratorBase
    from .md_system_base import MDStateDict
    
    logger = logging.getLogger(__name__)
    CreateDynFunc = Callable[[Atoms], Dynamics]
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator)class ASEIntegrator(MDIntegratorBase):
        """Base class for the integrators using ASE backend."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator.__init__)    def __init__(self, get_dynamics_fn: CreateDynFunc, timestep: float) -> None:
            """Initialize the integrator using ASE backend.
    
            Args:
                get_dynamics_fn (CreateDynFunc):
                  Function object for creating the Dynamics object of ASE.
                timestep (float): Integration time step in fs unit.
            """
            super().__init__(timestep)
            self._create_dyn = get_dynamics_fn
            self._dyn = None
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator.create_ase_dynamics)    def create_ase_dynamics(self, atoms: Atoms) -> Dynamics:
            """Create the ASE's Dynamics object.
    
            Args:
                atoms (Atoms): ASE's Atoms object containing the system to simulate.
            Returns:
                Dynamics : ASE's Dynamics object.
            """
            self._dyn = self._create_dyn(atoms)
            return self._dyn
    
    
    
        @property
        def ase_dynamics(self) -> Dynamics:
            """Get the ASE's Dynamics object.
    
            This method returns None, if ASE's dynamics object
            is not created by create_ase_dynamics() method
    
            Returns:
                Dynamics : ASE's Dynamics object.
            """
            return self._dyn
    
        @property
        def current_step(self) -> int:
            """Get the current integration step.
    
            Returns:
                int : Number of the current integration step.
            """
            result: int = self.ase_dynamics.get_number_of_steps()
            return result
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.ASEIntegrator.html#matlantis_features.features.md.ase_integrators.ASEIntegrator.set_temperature)    def set_temperature(self, value: float) -> None:
            """Set the target temperature of integrator.
    
            If this integrator does not support temperature control,
            this method raises the NotImplementedError exception.
    
            Args:
                value (float): Temperature in K unit.
            """
            self.ase_dynamics.set_temperature(temperature_K=value)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.html#matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator)class VelocityVerletIntegrator(ASEIntegrator):
        """Velocity Verlet integrator using ASE backend."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.html#matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.__init__)    def __init__(self, timestep: float) -> None:
            """Create a Velocity Verlet integrator.
    
            Args:
                timestep (float): Integration time step in fs unit.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                return VelocityVerlet(atoms, self.timestep * units.fs)
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.html#matlantis_features.features.md.ase_integrators.VelocityVerletIntegrator.set_temperature)    def set_temperature(self, value: float) -> None:
            """Set the target temperature of integrator.
    
            This integrator does not support temperature control,
            so this method always raises the NotImplementedError exception.
    
            Args:
                value (float): Temperature in K unit.
            """
            raise RuntimeError("this ensemble does not support temperature control")
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.LangevinIntegrator.html#matlantis_features.features.md.ase_integrators.LangevinIntegrator)class LangevinIntegrator(ASEIntegrator):
        """Langevin dynamics integrator generating NVT ensemble, using the ASE backend."""
    
        rng_state: Optional[np.random.RandomState]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.LangevinIntegrator.html#matlantis_features.features.md.ase_integrators.LangevinIntegrator.__init__)    def __init__(
            self,
            timestep: float,
            temperature: float,
            friction: float = 0.1,
            fixcm: bool = True,
            random_seed: int = 0,
        ) -> None:
            """Create Langevin dynamics integrator.
    
            Args:
                timestep (float): Integration time step in fs unit.
                temperature (float): Target temperature of the integrator in K unit.
                friction (float, optional):
                  Friction coefficient in the Langevin equation.
                  Defaults to 0.1.
                fixcm (bool, optional):
                  If True, the position and momentum of the center of mass is kept unperturbed.
                  Defaults to True.
                random_seed (int, optional):
                  Random number generator's seed used for the simulation.
                  Defaults to 0.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                if self.rng_state is None:
                    rng = np.random.RandomState(random_seed)
                else:
                    rng = self.rng_state
                    self.rng_state = None
                dyn = Langevin(
                    atoms,
                    timestep=self.timestep * units.fs,
                    temperature_K=self.temperature,
                    friction=self.friction,
                    fixcm=self.fixcm,
                    rng=rng,
                )
                return dyn
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
            self.temperature = temperature
            self.friction = friction
            self.fixcm = fixcm
            self.rng_state = None
    
    
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the integrator.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            dyn = self.ase_dynamics
            return {"rng": dyn.rng}
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the integrator.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            self.rng_state = state["rng"]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.AndersenIntegrator.html#matlantis_features.features.md.ase_integrators.AndersenIntegrator)class AndersenIntegrator(ASEIntegrator):
        """Andersen thermostat integrator generating NVT ensemble, using the ASE backend."""
    
        rng_state: Optional[np.random.RandomState]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.AndersenIntegrator.html#matlantis_features.features.md.ase_integrators.AndersenIntegrator.__init__)    def __init__(
            self,
            timestep: float,
            temperature: float,
            andersen_prob: float = 0.01,
            fixcm: bool = True,
            random_seed: int = 0,
        ):
            """Create the Andersen integrator.
    
            Args:
                timestep (float): Integration time step in fs unit.
                temperature (float): Target temperature of the integrator in K unit.
                andersen_prob (float, optional): A random collision probability (0-1).
                  with this probability, atoms get assigned random velocity components.
                  Defaults to 0.01.
                fixcm (bool, optional):
                  If True, the position and momentum of the center of mass is kept unperturbed.
                  Defaults to True.
                random_seed (int, optional):
                  Random number generator's seed used for the simulation.
                  Defaults to 0.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                if self.rng_state is None:
                    rng = np.random.RandomState(random_seed)
                else:
                    rng = self.rng_state
                    self.rng_state = None
                dyn = Andersen(
                    atoms,
                    timestep=self.timestep * units.fs,
                    temperature_K=self.temperature,
                    andersen_prob=self.andersen_prob,
                    fixcm=self.fixcm,
                    rng=rng,
                )
                return dyn
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
            self.temperature = temperature
            self.andersen_prob = andersen_prob
            self.fixcm = fixcm
            self.rng_state = None
    
    
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the integrator.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            dyn = self.ase_dynamics
            return {"rng": dyn.rng}
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the integrator.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            self.rng_state = state["rng"]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator.html#matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator)class NVTBerendsenIntegrator(ASEIntegrator):
        """Berendsen thermostat integrator generating NVT ensemble, using the ASE backend."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator.html#matlantis_features.features.md.ase_integrators.NVTBerendsenIntegrator.__init__)    def __init__(
            self,
            timestep: float,
            temperature: float,
            taut: float = 500,
            fixcm: bool = True,
        ):
            """Create the Berendsen thermostat integrator.
    
            Args:
                timestep (float): Integration time step in fs unit.
                temperature (float): Target temperature of the integrator in K unit.
                taut (float, optional):
                  Time constant for Berendsen temperature coupling in fs.
                  Defaults to 500 fs.
                fixcm (bool, optional):
                  If True, the position and momentum of the center of mass is kept unperturbed.
                  Defaults to True.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                return NVTBerendsen(
                    atoms,
                    timestep=self.timestep * units.fs,
                    temperature_K=self.temperature,
                    taut=self.taut * units.fs,
                    fixcm=self.fixcm,
                )
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
            self.temperature = temperature
            self.taut = taut
            self.fixcm = fixcm
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator.html#matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator)class NPTBerendsenIntegrator(ASEIntegrator):
        """Berendsen thermo/barostat integrator generating an NPT ensemble via ASE.
    
        Relax both temperature and isotropic pressure.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator.html#matlantis_features.features.md.ase_integrators.NPTBerendsenIntegrator.__init__)    def __init__(
            self,
            timestep: float,
            temperature: float,
            pressure: float,
            taut: float = 500,
            taup: float = 1000,
            compressibility: float = 67.2,
            fixcm: bool = True,
        ):
            """Create Berendsen thermo/barostat integrator.
    
            Args:
                timestep (float): Integration time step in fs.
                temperature (float): Target temperature in K.
                pressure (float): Target isotropic pressure in eV/Angstrom^3.
                taut (float, optional):
                  Temperature coupling time constant (fs). Smaller values enforce temperature more strongly; defaults to 500 fs.
                taup (float, optional):
                  Time constant for Berendsen pressure coupling.  Defaults to 1 ps.
                compressibility (float, optional):
                  The compressibility of the material (Angstrom^3/eV).
                  Defaults to 67.2.
                fixcm (bool, optional):
                  If True, the position and momentum of the center of mass is kept unperturbed.
                  Defaults to True.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                return NPTBerendsen(
                    atoms,
                    timestep=self.timestep * units.fs,
                    temperature_K=self.temperature,
                    pressure_au=self.pressure,
                    taut=self.taut * units.fs,
                    taup=self.taup * units.fs,
                    compressibility_au=self.compressibility,
                    fixcm=self.fixcm,
                )
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
            self.temperature = temperature
            self.pressure = pressure
            self.taut = taut
            self.taup = taup
            self.compressibility = compressibility
            self.fixcm = fixcm
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NPTIntegrator.html#matlantis_features.features.md.ase_integrators.NPTIntegrator)class NPTIntegrator(ASEIntegrator):
        """Nose-Hoover/Parrinello-Rahman dynamics integrator generating NPT ensemble,\
        using the ASE backend."""
    
        npt_state: Optional[Dict[str, Any]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.md.ase_integrators.NPTIntegrator.html#matlantis_features.features.md.ase_integrators.NPTIntegrator.__init__)    def __init__(
            self,
            timestep: float,
            temperature: float,
            pressure: float,
            ttime: float,
            pfactor: float,
            mask: Optional[np.ndarray] = None,
            hydrostatic_strain: bool = False,
        ):
            """Create the Nose-Hoover/Parrinello-Rahman dynamics integrator.
    
            Args:
                timestep (float): Integration time step in fs unit.
                temperature (float): Target temperature of the integrator in K unit.
                pressure (float): Target pressure  in eV/Angstrom^3 unit.
                ttime (float): Characteristic timescale of the thermostat, in fs units.
                pfactor (float): A constant in the barostat differential equation.
                  If a characteristic barostat timescale of ptime is desired,
                  set pfactor to ptime^2 * B (where ptime is in units matching eV, Angstrom, u;
                  and B is the Bulk Modulus, given in eV/Angstrom^3).
                mask (np.ndarray or None, optional):
                  A 3x3 array indicating if the system can change size along the three Cartesian axes.
                  Defaults to None.
                hydrostatic_strain (bool, optional):
                  Constrain the cell by only allowing hydrostatic deformation.
                  The virial tensor is replaced by np.diag([np.trace(virial)]*3). Defaults to False.
            """
    
            def _get_dyn(atoms: Atoms) -> Dynamics:
                if hydrostatic_strain:
                    atoms = IsotropicFilter(atoms)
                dyn = NPT(
                    atoms,
                    timestep=self.timestep * units.fs,
                    temperature_K=self.temperature,
                    externalstress=self.pressure,
                    ttime=self.ttime * units.fs,
                    pfactor=self.pfactor,
                    mask=self.mask,
                )
                if self.npt_state is not None:
                    for k, v in self.npt_state.items():
                        setattr(dyn, k, v)
                    self.npt_state = None
                return dyn
    
            super().__init__(get_dynamics_fn=_get_dyn, timestep=timestep)
            self.temperature = temperature
            self.pressure = pressure
            self.ttime = ttime
            self.pfactor = pfactor
            self.mask = mask
            self.hydrostatic_strain = hydrostatic_strain
            self.npt_state = None
    
    
    
        @property
        def state(self) -> MDStateDict:
            """Get the serializable state of the integrator.
    
            Returns:
                MDStateDict : Serializable state data.
            """
            dyn = self.ase_dynamics
            state = dyn.get_data()
            # print(state)
            return {"npt_state": state}
    
        @state.setter
        def state(self, state: MDStateDict) -> None:
            """Set the serializable state of the integrator.
    
            Args:
                state (MDStateDict): Serializable state data.
            """
            self.npt_state = state["npt_state"]
    
    
    
