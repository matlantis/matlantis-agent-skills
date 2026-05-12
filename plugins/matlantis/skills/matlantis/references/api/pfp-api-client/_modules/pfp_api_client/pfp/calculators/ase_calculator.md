# Source code for pfp_api_client.pfp.calculators.ase_calculator
    
    
    import threading
    from types import TracebackType
    from typing import Dict, List, Optional, Sequence, Set, Type, cast
    
    import numpy as np
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes
    
    from pfp_api_client import Estimator
    from pfp_api_client.logging import get_logger
    from pfp_api_client.pfp.estimator import EstimatorSystem
    from pfp_api_client.pfp.utils.exceptions import ConcurrentUseDetected
    from pfp_api_client.pfp.utils.interrupt_watcher import InterruptWatcher
    from pfp_api_client.pfp.utils.messages import MessageEnum, MessageType, message_info
    
    logger = get_logger(__name__)
    
    _logging_type = {
        MessageType.CustomMessage: logger.info,
        MessageType.CustomWarning: logger.warning,
    }
    
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.calculators.html#pfp_api_client.pfp.calculators.ase_calculator.ASECalculator)class ASECalculator(Calculator):  # type: ignore
        """
        ASECalculator is a derivation of `ase.calculators.calculator`.
        Usually this will be used from `ase.Atoms` object.
        It can be safe to have multiple instances of `ASECalculator` running
        concurrently in different threads, as long as they do not share
        the same `Estimator` instance. It is *not* safe to use a single
        instance of `ASECalculator` concurrently and trying to do so will
        likely raise a `ConcurrentUseDetected` exception.
    
        Parameters
        --------
        estimator: Estimator
            PFP Estimator
        """
    
        atoms: Optional[Atoms]
        estimator: Estimator
        implemented_properties: List[str]
        current_messages: Set[MessageEnum]
    
        def __init__(self, estimator: Estimator) -> None:
            super(ASECalculator, self).__init__()
            estimator.attach_to_calculator(self)  # required to avoid accidental one-to-many usage
            self.estimator = estimator
            self.implemented_properties = self.convert_properties(estimator.implemented_properties)
            self.current_messages = set()
            self.default_properties: List[str] = [
                prop
                for prop in ["energy", "charges", "forces", "free_energy"]
                if prop in self.implemented_properties
            ]
            self.concurrency_checker = threading.Lock()
            self._last_descriptor_indices: Optional[np.ndarray] = None
    
        def __enter__(self) -> "ASECalculator":
            return self
    
        def __exit__(
            self,
            exc_type: Optional[Type[BaseException]],
            exc_val: Optional[BaseException],
            exc_tb: Optional[TracebackType],
        ) -> None:
            self.estimator.close()
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.calculators.html#pfp_api_client.pfp.calculators.ase_calculator.ASECalculator.set_default_properties)    def set_default_properties(self, default_properties: Sequence[str]) -> None:
            """
            Specify the properties that are always calculated.
    
            Parameters
            --------
            default_properties: Sequence[str]
                Properties to always calculate. Must be a subset of ``implemented_properties``.
            """
            for prop in default_properties:
                assert prop in self.implemented_properties
    
            # make a copy to prevent issues with outside modifications of the parameter
            self.default_properties = list(default_properties)
    
    
    
        def convert_properties(self, properties: Sequence[str]) -> List[str]:
            properties_for_calc = []
            for property in properties:
                if property == "virial":
                    properties_for_calc.append("stress")
                elif property == "energy":
                    properties_for_calc.append(property)
                    # NOTE: `FrechetCellFilter` , `ExpCellFilter` and `UnitCellFilter` require
                    # a calculator to be able to calculate free_energy property
                    properties_for_calc.append("free_energy")
                else:
                    properties_for_calc.append(property)
            return properties_for_calc
    
        def reverse_convert_properties(self, properties: Sequence[str]) -> List[str]:
            def convert(s: str) -> str:
                if s == "stress":
                    return "virial"
                else:
                    return s
    
            return [convert(s) for s in properties if s != "free_energy"]
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.calculators.html#pfp_api_client.pfp.calculators.ase_calculator.ASECalculator.calculate)    def calculate(
            self,
            atoms: Optional[Atoms] = None,
            properties: Optional[Sequence[str]] = None,
            system_changes: Sequence[str] = all_changes,
        ) -> None:
            """
            Do the calculation. This function will be called from various ASE objects such as \
            `ase.optimizer` or `ase.md` classes.
    
            Parameters
            --------
            atoms:
                ase.Atoms object.
            properties:
                Sequence of what needs to be calculated. Current supported values are 'energy', \
                'forces', 'stress', and 'charges'.
            system_changes:
                Sequence of what has changed since last calculation.
            """
    
            # try to get a lock, and if it fails then make it clear that this
            # usage is forbidden, we must not use one ASECalculator concurrently
            got_lock = self.concurrency_checker.acquire(blocking=False)
            if not got_lock:
                raise ConcurrentUseDetected()
            try:
                # check that we were not interrupted by SIGINT
                InterruptWatcher.assert_not_interrupted()
                return self._calculate(atoms, properties, system_changes)
            finally:
                self.concurrency_checker.release()
    
    
    
        def _calculate(
            self,
            atoms: Optional[Atoms],
            properties: Optional[Sequence[str]],
            system_changes: Sequence[str],
        ) -> None:
            if properties is None:
                properties = ["energy"]
    
            # NOTE for concurrency: this method *must not* modify any of the
            #      parameters (we don't "own" them) and ideally should work
            #      on copies if we want to ensure consistency during the lifetime
            #      of this function.
    
            # the superclass only creates a copy of `atoms` and doesn't touch the other parameters
            super(ASECalculator, self).calculate(atoms, properties, system_changes)
            assert self.atoms is not None
            assert id(self.atoms) != id(atoms)  # superclass should make a copy
    
            pbc = self.atoms.get_pbc().astype(np.uint8)
    
            _properties = list(properties)  # work on a copy of the passed value
            if len(self.default_properties) > 0:
                for prop in self.default_properties:
                    if prop not in _properties:
                        _properties.append(prop)
    
            if (
                "forces" in _properties
                and np.all(pbc)
                and "stress" not in _properties
                and "stress" in self.implemented_properties
            ):
                # Automatically add stress when force is requested
                _properties.append("stress")
            if self.atoms.cell.rank < np.sum(pbc):
                logger.warning(
                    "Cell seems to be not defined since %d lattice vectors and its volume not defined",
                    self.atoms.cell.rank,
                )
            cell = self.atoms.get_cell(complete=True)
            cell_array = cell.array
    
            atomic_numbers = self.atoms.get_atomic_numbers()
            coordinates = self.atoms.get_positions()
            properties_estimator = self.reverse_convert_properties(_properties)
    
            self.results = self.estimator.estimate(
                EstimatorSystem(
                    properties=properties_estimator,
                    coordinates=coordinates,
                    cell=cell_array,
                    atomic_numbers=atomic_numbers,
                    pbc=pbc,
                )
            )
    
            # NOTE (himkt): to make compatible with ``SinglePointCalculator``.
            for key, value in self.results.items():
                if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.floating):
                    self.results[key] = value.astype(np.float64)
    
            if "energy" in self.results:
                self.results["free_energy"] = self.results["energy"]
    
            if "virial" in self.results and np.all(pbc):
                if not isinstance(self.results["virial"], np.ndarray):
                    raise ValueError(
                        f"received virial must be np.ndarray, but got {type(self.results['messages'])}"
                    )
                self.results["stress"] = self.calculate_stress(
                    self.results["virial"],
                    cell.volume,
                )
                del self.results["virial"]
    
            if "messages" in self.results:
                if not isinstance(self.results["messages"], list):
                    raise ValueError(
                        f"received messages must be list, but got {type(self.results['messages'])}"
                    )
                for message in self.results["messages"]:
                    if self.estimator.message_isactive[message]:
                        _logging_type[message_info(message)[0]](str(message))
                        self.current_messages.add(message)
                del self.results["messages"]
    
        @staticmethod
        def calculate_stress(virial: np.ndarray, atom_volume: np.ndarray) -> np.ndarray:
            """Calculate stress from virial and cell volume."""
            stress: np.ndarray = virial / atom_volume
            return stress
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.calculators.html#pfp_api_client.pfp.calculators.ase_calculator.ASECalculator.get_descriptors)    def get_descriptors(
            self, atoms: Atoms, indices: Optional[np.ndarray] = None
        ) -> Dict[str, np.ndarray]:
            """
            Calculate and return scalar descriptors for a given set of atoms.
    
            Parameters
            ----------
            atoms : ase.Atoms
                The ASE Atoms object for which descriptors are to be calculated.
            indices : Optional[np.ndarray], optional
                A 1-dimensional array of indices (between 0 and n_atoms-1) specifying the
                atomic indices for descriptors to be obtained.
    
            Returns
            -------
            descriptors : Dict[str, np.ndarray]
                Dictionary containing the calculated scalar descriptors under the key `scalar`.
                The descriptors have a shape of [n_atoms, 256].
            """
    
            def _indices_changed(
                current: Optional[np.ndarray], previous: Optional[np.ndarray]
            ) -> bool:
                # Both None: no change
                if current is None and previous is None:
                    return False
                # One is None, the other is not: changed
                elif current is None or previous is None:
                    return True
                # Both are arrays: compare contents
                else:
                    return not np.array_equal(current, previous)
    
            if not self.calculation_required(atoms, ["descriptors"]) and not _indices_changed(
                indices, self._last_descriptor_indices
            ):
                return cast(Dict[str, np.ndarray], self.results["descriptors"])
    
            self._last_descriptor_indices = indices.copy() if indices is not None else None
    
            # NOTE (mao): self.atoms is used for comparison in the `calculation_required` method
            #             defined in the `Calculator` superclass.
            self.atoms = atoms.copy()
    
            pbc = atoms.get_pbc().astype(np.uint8)
            cell = atoms.get_cell(complete=True)
            cell_array = cell.array
            atomic_numbers = atoms.get_atomic_numbers()
            coordinates = atoms.get_positions()
    
            results = self.estimator.get_descriptors(
                EstimatorSystem(
                    properties=["energy"],
                    coordinates=coordinates,
                    cell=cell_array,
                    atomic_numbers=atomic_numbers,
                    pbc=pbc,
                ),
                indices=indices,
            )
            # NOTE (himkt): to make compatible with ``SinglePointCalculator``.
            for key, value in results.items():
                if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.floating):
                    results[key] = value.astype(np.float64)
    
            if "messages" in results:
                if not isinstance(results["messages"], list):
                    raise ValueError(
                        f"received messages must be list, but got {type(results['messages'])}"
                    )
                for message in results["messages"]:
                    if self.estimator.message_isactive[message]:
                        _logging_type[message_info(message)[0]](str(message))
                        self.current_messages.add(message)
                del results["messages"]
    
            del results["calc_stats"]
    
            descriptors: Dict[str, np.ndarray] = {"scalar": cast(np.ndarray, results["descriptors"])}
            self.results["descriptors"] = descriptors  # type: ignore
    
            return descriptors
    
    
    
