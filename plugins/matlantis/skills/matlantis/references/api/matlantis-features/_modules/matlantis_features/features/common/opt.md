# Source code for matlantis_features.features.common.opt
    
    
    import contextlib
    import logging
    from copy import copy, deepcopy
    from dataclasses import dataclass
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union
    
    import ase
    import numpy as np
    from packaging import version
    
    if version.parse(ase.__version__) >= version.parse("3.23.0"):
        from ase.filters import Filter
        from ase.mep import NEB
    else:
        from ase.constraints import Filter
        from ase.neb import NEB
    
    from ase.optimize import BFGS, FIRE, LBFGS, BFGSLineSearch, LBFGSLineSearch, MDMin
    from ase.optimize.optimize import Optimizer
    from tqdm.auto import tqdm
    from tqdm.contrib.logging import logging_redirect_tqdm
    
    from matlantis_features.ase_ext.optimize import FIRELBFGS
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase, FeatureBaseOutput, FeatureBaseResult
    from matlantis_features.filters import MatlantisFilter, UnitCellASEFilter
    from matlantis_features.utils import FeatureCost, get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    from matlantis_features.utils.save_util import (
        get_calculator_info,
        get_user_info,
        get_version_info,
        serialize_single_item,
    )
    
    _LATER_OR_EQUALS_TO_ASE_3_26_0 = version.parse(ase.__version__) >= version.parse("3.26.0")
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.OptFeatureOutput.html#matlantis_features.features.common.opt.OptFeatureOutput)@dataclass
    class OptFeatureOutput(FeatureBaseOutput):
        """Class containing the outputs for results of various optimizer features."""
    
        # TODO: replace MatlantisAtoms
        atoms: MatlantisAtoms
        converged: bool
        steps: int
        # TODO: add force_max
        # force_max: float
        energy_log: List[float]
        force_log: List[float]
        description: str
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.OptFeatureResult.html#matlantis_features.features.common.opt.OptFeatureResult)@dataclass
    class OptFeatureResult(FeatureBaseResult):
        """Class with the results for various optimizer features."""
    
        feature: "OptFeatureBase"
        output: "OptFeatureOutput"
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase)class OptFeatureBase(FeatureBase):
        """Template base class for optimizer features.
    
        This class is inherited by other optimizer features, \
        and not intended to be used directly.
        """
    
        estimator_fn: Optional[EstimatorFnType]
    
        def __call__(self, atoms: Union[ASEAtoms, MatlantisAtoms]) -> OptFeatureResult:
            """Call function.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms): The system to be optimized.
            Returns:
                OptFeatureResult :
                    A dataclass containing information about the result of optimization.
                    Raises a NotImplementedError in the base class.
            """
            raise NotImplementedError
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for Used for keeping track of cost information.
    
            Args:
                n_atoms (int): The number of atoms.
            Returns:
                FeatureCost :
                    A dataclass containing cost information for the feature.
                    Raises a NotImplementedError in the base class.
            """
            raise NotImplementedError
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the OptFeatureBase.
    
            Returns:
                dict[str, Any] : A dict containing a serialized OptFeatureBase.
            """
            raise NotImplementedError
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature)class ASEOptFeature(OptFeatureBase):
        """The matlantis-feature that wraps around an ASE Optimizer object."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature.__init__)    def __init__(
            self,
            get_optimizer: Callable[[ASEAtoms], Optimizer],
            n_run: int = 200,
            fmax: float = 0.001,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            estimator_fn: Optional[EstimatorFnType] = None,
        ):
            """Initialize an instance.
    
            Args:
                get_optimizer (ASEAtoms -> Optimizer): A function to get optimizer.
                n_run (int, optional):
                    The maximum number of optimization steps performed to try to
                    reach force convergence (specified by fmax). Defaults to 200.
                fmax (float, optional):
                    The maximum force (in eV/Angstrom) for convergence of optimization.
                    Defaults to 0.01.
                attach_methods (list[tuple[ -> None, int]] or None, optional):
                    Functions to be called by the optimizer. A list of tuples (or None) must be specified.
                    The first element of each tuple is a callable function, and the second element
                    is an integer that specifies the interval (in optimization steps) of when
                    the optimizer will call the function. If None, no functions will be attached.
                    Defaults to None.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 10.
                filter (bool or MatlantisFilter):
                    The MatlantisFilter used in the optimization. If True, the default filter, i.e.
                    UnitCellASEFilter, will be used for the cell shape optimization. If False, no
                    filter is used.
                    Note: The filter is only applicable in the optimization of ASEAtoms and MatlantisAtoms,
                    not NEB class.
                    Defaults to False.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            super(ASEOptFeature, self).__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            with self.init_scope():
                self.get_optimizer = get_optimizer
                self.n_run = n_run
                self.fmax = fmax
                self.attach_methods: List[Tuple[Callable[[], None], int]]
                if attach_methods:
                    self.attach_methods = attach_methods
                else:
                    self.attach_methods = []
    
                self.show_progress_bar = show_progress_bar
                self.tqdm_options = tqdm_options
                self.show_logger = show_logger
                self.logger_interval = logger_interval
                self.filter: Optional[MatlantisFilter] = None
                if filter is True:
                    self.filter = UnitCellASEFilter()
                elif isinstance(filter, MatlantisFilter):
                    self.filter = filter
                self.estimator_fn = estimator_fn
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature.attach)    def attach(self, extn: Callable[[], None], interval: int) -> None:
            """Used for attaching a function to the ASE optimizer.
    
            Args:
                extn ( -> None): The function to be attached.
                interval (int):
                    The interval (in optimization steps) of when the optimizer will call the function.
            """
            self.attach_methods.append((extn, interval))
    
    
    
        def _check_input(self, system: Union[ASEAtoms, NEB, Filter]) -> None:
            if isinstance(system, Filter):
                logging.warning(
                    "The support of ASE Filter class will be deprecated in the future. "
                    "Please use the 'filter' parameter in the __init__ function."
                )
            elif isinstance(system, NEB) and self.filter:
                logging.warning("Cannot optimize NEB together with filter. The setting of 'filter' will be ignored.")
    
        def __call__(self, atoms: Union[ASEAtoms, MatlantisAtoms, NEB, Filter]) -> OptFeatureResult:
            """Call function.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms or NEB or Filter):
                    The system to be optimized. In addition to ASEAtoms and MatlantisAtoms,
                    NEB (Matlantis) and Filter (ASE) objects are also supported.
                    Note: the support of ASE Filter class will be deprecated in the future.
                    Please use the 'filter' parameter in the __init__ function.
            Returns:
                OptFeatureResult : A dataclass containing information about the result of optimization.
            """
            logger = logging.getLogger(__name__)
            system: ASEAtoms
            if isinstance(atoms, MatlantisAtoms):
                system = atoms.ase_atoms
            else:
                system = atoms
    
            self._check_input(system)
    
            with contextlib.ExitStack() as stack:
                # TODO(tgp) allow the user to set a calculator (i.e., don't overwrite it)
                filter_wrapped = False
                if isinstance(system, ASEAtoms):  # type: ignore
                    system.calc = stack.enter_context(get_calculator(estimator_fn=self.estimator_fn))
                    if self.filter is not None:
                        system = self.filter(system)
                        filter_wrapped = True
                elif isinstance(system, Filter):  # deprecated
                    system.atoms.calc = stack.enter_context(get_calculator(estimator_fn=self.estimator_fn))
                elif isinstance(system, NEB):
                    if len(system.images) != len(set(id(img) for img in system.images)):
                        raise ValueError("Identical Atoms instances are not allowed in the NEB input.")
    
                    # always attach a different Calculator to each image even if parallel=False,
                    # because otherwise frequent cache invalidation can make computation very slow
                    for image in system.images:
                        image.calc = stack.enter_context(get_calculator(estimator_fn=self.estimator_fn))
    
                    if not system.parallel:
                        logger.info("NEB computation is NOT using parallelization.")
    
                opt = self.get_optimizer(system)
    
                energy_log: List[float] = list()
                force_log: List[float] = list()
    
                pbar = None
                if self.show_progress_bar:
                    ntotal = self.n_run + 1
                    if self.tqdm_options is None:
                        pbar = tqdm(total=ntotal)
                    else:
                        pbar = tqdm(total=ntotal, **self.tqdm_options)
    
                def update_log() -> None:
                    E = system.get_potential_energy()
                    F = np.max(np.linalg.norm(system.get_forces(), axis=1))
                    energy_log.append(float(E))
                    force_log.append(float(F))
                    with logging_redirect_tqdm():
                        if self.show_logger and (opt.nsteps % self.logger_interval == 0):
                            if isinstance(system, NEB):
                                energies: List[float] = [float(atoms.get_potential_energy()) for atoms in system.images]
                                logger.info(
                                    f"steps: {opt.nsteps:5d}  "
                                    f"energies: {' '.join([f'{e:5.2f}' for e in energies])} eV  "
                                    f"max_force: {F:5.2f} eV/Ang"
                                )
                            else:
                                logger.info(f"steps: {opt.nsteps:5d}  energy：{E / len(system):5.2f} eV/atom  max_force: {F:5.2f} eV/Ang")
                        if pbar is not None:
                            pbar.update(1)
    
                opt.attach(update_log, interval=1)
                # TODO: add interval info?
                if self.attach_methods:
                    for extn, interval in self.attach_methods:
                        opt.attach(extn, interval=interval)
                opt.run(steps=self.n_run, fmax=self.fmax)
                converged = _get_convergence_status(opt)
                if pbar is not None:
                    pbar.update(self.n_run - opt.nsteps)
                    pbar.close()
    
            if filter_wrapped:
                system = system.atoms
    
            output = OptFeatureOutput(
                atoms=MatlantisAtoms(system),
                converged=converged,
                steps=opt.get_number_of_steps(),
                # TODO: add force_max
                # "force_max": force_log[-1],
                energy_log=energy_log,
                force_log=force_log,
                description="",
            )
    
            return OptFeatureResult(
                feature=deepcopy(self),
                call_params={"atoms": atoms.copy() if isinstance(atoms, MatlantisAtoms) else copy(atoms)},
                output=output,
                calculator_info=get_calculator_info(estimator_fn=self.estimator_fn),
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
            return FeatureCost(n_pfp_run=self.n_run, time_additional=0.1)
    
    
    
    
    _OPT_DOC_STRING_WITHOUT_MAX_STEP = """Initialize an instance.
    
            Args:
                n_run (int, optional):
                    The maximum number of optimization steps performed to try to
                    reach force convergence (specified by fmax). Defaults to 200.
                fmax (float, optional):
                    The maximum force (in eV/Angstrom) for convergence of optimization.
                    Defaults to 0.01.
                attach_methods (list[tuple[ -> None, int]] or None, optional):
                    Functions to be called by the optimizer. A list of tuples must be specified,
                    where the first element of each tuple is a callable function, and the second element
                    is an integer N. The optimizer will call the function every N steps.
                    If None, no functions will be attached. Defaults to None.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 10.
                filter (bool or MatlantisFilter, optional):
                    The MatlantisFilter used in the optimization. If True, the default filter, i.e.
                    UnitCellASEFilter, will be used for the cell shape optimization. If False, no
                    filter is used.
                    Note: The filter is only applicable in the optimization of ASEAtoms and MatlantisAtoms,
                    not NEB class.
                    Defaults to False.
                trajectory (str or None, optional):
                    Filename for trajectory file. If None, no trajectory file will be saved. Defaults
                    to None.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
    _OPT_DOC_STRING = (
        _OPT_DOC_STRING_WITHOUT_MAX_STEP
        + """
                maxstep (float, optional):
                    The maximum distance an atom can move per iteration. The unit is angstrom.
                    Defaults to 0.2.
            """
    )
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.MDMinASEOptFeature.html#matlantis_features.features.common.opt.MDMinASEOptFeature)class MDMinASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the ASE implementation of the MDMin algorithm."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.MDMinASEOptFeature.html#matlantis_features.features.common.opt.MDMinASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
        ) -> None:
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return MDMin(atoms, trajectory=trajectory, logfile=None)
    
            super(MDMinASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING_WITHOUT_MAX_STEP
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.MDMinASEOptFeature.html#matlantis_features.features.common.opt.MDMinASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the MDMinASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized MDMinASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireASEOptFeature.html#matlantis_features.features.common.opt.FireASEOptFeature)class FireASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the Fast Inertial Relaxation Engine (FIRE), \
        as implemented within ASE.
    
        For details of the FIRE algorithm, please see:
        Guénolé et al., Comput. Mater. Sci. 175 (2020) 109584.
        """
    
        # TODO: add Fire related arguments
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireASEOptFeature.html#matlantis_features.features.common.opt.FireASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep: float = 0.2,
        ):
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
            self.maxstep = maxstep
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return FIRE(atoms, maxstep=maxstep, trajectory=trajectory, logfile=None)
    
            super(FireASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireASEOptFeature.html#matlantis_features.features.common.opt.FireASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the FireASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized FireASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep": self.maxstep,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSASEOptFeature.html#matlantis_features.features.common.opt.BFGSASEOptFeature)class BFGSASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the ASE implementation of the BFGS algorithm."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSASEOptFeature.html#matlantis_features.features.common.opt.BFGSASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep: float = 0.2,
        ):
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
            self.maxstep = maxstep
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return BFGS(atoms, maxstep=maxstep, trajectory=trajectory, logfile=None)
    
            super(BFGSASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSASEOptFeature.html#matlantis_features.features.common.opt.BFGSASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the BFGSASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized BFGSASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep": self.maxstep,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSASEOptFeature.html#matlantis_features.features.common.opt.LBFGSASEOptFeature)class LBFGSASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the ASE implementation of the LBFGS algorithm."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSASEOptFeature.html#matlantis_features.features.common.opt.LBFGSASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep: float = 0.2,
        ):
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
            self.maxstep = maxstep
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return LBFGS(atoms, maxstep=maxstep, trajectory=trajectory, logfile=None)
    
            super(LBFGSASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSASEOptFeature.html#matlantis_features.features.common.opt.LBFGSASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the LBFGSASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized LBFGSASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep": self.maxstep,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature)class BFGSLineSearchASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the ASE implementation of the BFGS line search algorithm."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep: float = 0.2,
        ):
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
            self.maxstep = maxstep
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return BFGSLineSearch(atoms, maxstep=maxstep, trajectory=trajectory, logfile=None)
    
            super(BFGSLineSearchASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the BFGSLineSearchASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized BFGSLineSearchASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep": self.maxstep,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature)class LBFGSLineSearchASEOptFeature(ASEOptFeature):
        """The matlantis-feature that uses the ASE implementation of the BFGS line search algorithm."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep: float = 0.2,
        ):
            self.n_run = n_run
            self.fmax = fmax
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            self.estimator_fn = estimator_fn
            self.maxstep = maxstep
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return LBFGSLineSearch(atoms, maxstep=maxstep, trajectory=trajectory, logfile=None)
    
            super(LBFGSLineSearchASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the LBFGSLineSearchASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized LBFGSLineSearchASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep": self.maxstep,
            }
            return res
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireLBFGSASEOptFeature.html#matlantis_features.features.common.opt.FireLBFGSASEOptFeature)class FireLBFGSASEOptFeature(ASEOptFeature):
        """The matlantis-feature that combines FIRE and LBFGS optimization algorithms."""
    
        # TODO: add Fire related arguments
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireLBFGSASEOptFeature.html#matlantis_features.features.common.opt.FireLBFGSASEOptFeature.__init__)    def __init__(
            self,
            n_run: int = 200,
            fmax: float = 0.01,
            switch: float = 0.05,
            switch_decrease_rate: float = 0.9,
            switch_patience_fire: int = 10,
            switch_patience_lbfgs: int = 10,
            attach_methods: Optional[List[Tuple[Callable[[], None], int]]] = None,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            logger_interval: int = 10,
            filter: Union[bool, MatlantisFilter] = False,
            trajectory: Optional[str] = None,
            estimator_fn: Optional[EstimatorFnType] = None,
            maxstep_fire: float = 0.2,
            maxstep_lbfgs: float = 0.2,
        ):
            """Initialize an instance.
    
            Args:
                n_run (int, optional):
                    The maximum number of optimization steps performed to try to
                    reach force convergence (specified by fmax). Defaults to 200.
                fmax (float, optional):
                    The maximum force (in eV/Angstrom) for convergence of optimization.
                    Defaults to 0.01.
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
                attach_methods (list[tuple[ -> None, int]] or None, optional):
                    Functions to be called by the optimizer. A list of tuples must be specified,
                    where the first element of each tuple is a callable function, and the second element
                    is an integer N. The optimizer will call the function every N steps.
                    If None, no functions will be attached. Defaults to None.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                logger_interval (int, optional):
                    The interval of when to print out logger information.
                    Defaults to 10.
                filter (bool or MatlantisFilter):
                    The MatlantisFilter used in the optimization. If True, the default filter, i.e.
                    UnitCellASEFilter, will be used for the cell shape optimization. If False, no
                    filter is used.
                    Note: The filter is only applicable in the optimization of ASEAtoms and MatlantisAtoms,
                    not NEB class.
                    Defaults to False.
                trajectory (str or None, optional):
                    Filename for trajectory file. If None, no trajectory file will be saved. Defaults
                    to None.
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
                maxstep_fire (float, optional):
                    The maximum distance an atom can move per iteration when using FIRE.
                    Defaults to 0.2.
                maxstep_lbfgs (float, optional):
                    The maximum distance an atom can move per iteration when using LBFGS.
                    Defaults to 0.2.
            """
            self.n_run = n_run
            self.fmax = fmax
            self.switch = switch
            self.switch_decrease_rate = switch_decrease_rate
            self.switch_patience_fire = switch_patience_fire
            self.switch_patience_lbfgs = switch_patience_lbfgs
            self.show_progress_bar = show_progress_bar
            self.tqdm_options = tqdm_options
            self.show_logger = show_logger
            self.logger_interval = logger_interval
            self.trajectory = trajectory
            # TODO: estimator_fn
            self.maxstep_fire = maxstep_fire
            self.maxstep_lbfgs = maxstep_lbfgs
    
            def get_optimizer(atoms: ASEAtoms) -> Optimizer:
                return FIRELBFGS(
                    atoms,
                    switch=switch,
                    switch_decrease_rate=switch_decrease_rate,
                    switch_patience_fire=switch_patience_fire,
                    switch_patience_lbfgs=switch_patience_lbfgs,
                    trajectory=trajectory,
                    maxstep_fire=maxstep_fire,
                    maxstep_lbfgs=maxstep_lbfgs,
                    logfile=None,
                )
    
            super(FireLBFGSASEOptFeature, self).__init__(
                get_optimizer=get_optimizer,
                n_run=n_run,
                fmax=fmax,
                attach_methods=attach_methods,
                show_progress_bar=show_progress_bar,
                tqdm_options=tqdm_options,
                show_logger=show_logger,
                logger_interval=logger_interval,
                filter=filter,
                estimator_fn=estimator_fn,
            )
    
    
    
        __init__.__doc__ = _OPT_DOC_STRING
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.common.opt.FireLBFGSASEOptFeature.html#matlantis_features.features.common.opt.FireLBFGSASEOptFeature.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the FireLBFGSASEOptFeature.
    
            Returns:
                dict[str, Any] : A dict containing a serialized FireLBFGSASEOptFeature.
            """
            res = {
                "cls_path": f"{self.__class__.__module__}.{self.__class__.__name__}",
                "n_run": self.n_run,
                "fmax": self.fmax,
                "switch": self.switch,
                "switch_decrease_rate": self.switch_decrease_rate,
                "switch_patience_fire": self.switch_patience_fire,
                "switch_patience_lbfgs": self.switch_patience_lbfgs,
                # attach_methods cannot be saved
                "show_progress_bar": self.show_progress_bar,
                "tqdm_options": self.tqdm_options,
                "show_logger": self.show_logger,
                "logger_interval": self.logger_interval,
                "filter": serialize_single_item(self.filter),
                "trajectory": self.trajectory,
                # TODO: estimator_fn cannot be saved
                "maxstep_fire": self.maxstep_fire,
                "maxstep_lbfgs": self.maxstep_lbfgs,
            }
    
            # TODO: implement helper function
            return res
    
    
    
    
    def _get_convergence_status(opt: Optimizer) -> bool:
        if _LATER_OR_EQUALS_TO_ASE_3_26_0:
            # ASE 3.26.0 changed the behavior of Optimizer.converged()
            gradient = opt.optimizable.get_gradient()
            return bool(opt.converged(gradient=gradient))
        else:
            return bool(opt.converged())
    
