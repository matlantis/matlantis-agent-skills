# Source code for pfp_api_client.pfp.estimator
    
    
    import os
    import threading
    import time
    from dataclasses import dataclass
    from enum import Enum
    from types import TracebackType
    from typing import (
        TYPE_CHECKING,
        Dict,
        List,
        Mapping,
        Optional,
        Sequence,
        Set,
        Tuple,
        Type,
        Union,
        cast,
    )
    
    import grpc
    import numpy as np
    from ase.calculators.calculator import Calculator
    from packaging import version
    
    from pfp_api_client import __version__
    from pfp_api_client._proto import pfp_pb2, pfp_pb2_grpc
    from pfp_api_client._proto.pfp_pb2 import (
        CALC_MODE_CRYSTAL,
        CALC_MODE_CRYSTAL_PLUS_D3,
        CALC_MODE_CRYSTAL_U0,
        CALC_MODE_CRYSTAL_U0_PLUS_D3,
        CALC_MODE_MOLECULE,
        CALC_MODE_PBE,
        CALC_MODE_PBE_PLUS_D3,
        CALC_MODE_PBE_U,
        CALC_MODE_PBE_U_PLUS_D3,
        CALC_MODE_R2SCAN,
        CALC_MODE_R2SCAN_PLUS_D3,
        CALC_MODE_UNKNOWN,
        CALC_MODE_WB97XD,
        ELEMENT_STATUS_EXPECTED,
        ELEMENT_STATUS_EXPERIMENTAL,
        METHOD_TYPE_AUTO,
        METHOD_TYPE_MNCORE,
        METHOD_TYPE_PFVM,
        METHOD_TYPE_PFVM_D3_PFVM,
        NdArray,
    )
    from pfp_api_client.logging import get_logger
    from pfp_api_client.pfp.utils.grpc_channel import ChannelSingleton
    from pfp_api_client.pfp.utils.messages import MessageEnum
    from pfp_api_client.pfp.utils.shift_energies import (
        gaussian_shift_energies,
        r2scan_wo_U_shift_energies_v1_8_0,
        r2scan_wo_U_shift_energies_v1_9_3,
        vasp_U_shift_energies_v1_0_0,
        vasp_U_shift_energies_v1_4_0,
        vasp_U_shift_energies_v1_7_1,
        vasp_wo_U_shift_energies_v1_0_0,
        vasp_wo_U_shift_energies_v1_4_0,
        vasp_wo_U_shift_energies_v1_7_1,
    )
    from pfp_api_client.utils.error_handler_client_interceptor import ErrorHandlerClientInterceptor
    from pfp_api_client.utils.messages import from_messages_pb
    from pfp_api_client.utils.ndarray import from_ndarray_pb, to_ndarray_pb
    from pfp_api_client.utils.retry_client_interceptor import RetryClientInterceptor
    
    from .config import Config
    from .utils.exceptions import MultiCalculatorUseDetected
    
    if TYPE_CHECKING:
        from pfp_api_client._proto.pfp_pb2 import CalcMode, MethodType, TargetElementStatus
    
    logger = get_logger(__name__)
    
    METADATA_REQUEST_ID = "x-request-id"
    METADATA_N_ATOMS = "n_atoms"
    METADATA_N_INDICES = "n_indices"
    METADATA_MODEL_VERSION = "model_version"
    METADATA_METHOD_TYPE = "method_type"
    METADATA_CALC_MODE = "calc_mode"
    METADATA_CLIENT_VERSION = "client_version"
    METADATA_CLIENT_SIDE_PRIORITY = "client-side-priority"
    METADATA_EXECUTION_CONTEXT = "execution-context"
    METADATA_APPLICATION_CONTEXT = "application-context"
    latest_version = Config.latest_model_version
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorCalcMode)class EstimatorCalcMode(Enum):
        """
        Enum class which is used to determine `calc_mode` in `Estimator`.
    
        :ivar R2SCAN: This mode corresponds to DFT calculations with r2SCAN functional and plane wave basis sets.
        :ivar R2SCAN_PLUS_D3: R2SCAN mode with D3 correction enabled. \
        This mode corresponds to DFT calculations with r2SCAN functional and plane wave basis sets, \
        plus a D3 dispersion correction.
        :ivar PBE: This mode corresponds to DFT calculations with PBE functional and plane wave basis sets.
        :ivar PBE_PLUS_D3: PBE mode with D3 correction enabled. \
        This mode corresponds to DFT calculations with plane wave basis sets, \
        plus a D3 dispersion correction.
        :ivar PBE_U: PBE mode with Hubbard U parameter. This mode corresponds \
         to DFT calculations with plane wave basis sets.
        :ivar PBE_U_PLUS_D3: PBE_U mode with D3 correction enabled. \
        This mode corresponds to DFT calculations with plane wave basis sets, \
        plus a D3 dispersion correction.
        :ivar WB97XD: This mode corresponds to DFT calculations with local basis \
        sets.
        """
    
        CRYSTAL = "crystal"
        CRYSTAL_PLUS_D3 = "crystal_plus_d3"
        CRYSTAL_U0 = "crystal_u0"
        CRYSTAL_U0_PLUS_D3 = "crystal_u0_plus_d3"
        MOLECULE = "molecule"
        R2SCAN = "r2scan"
        R2SCAN_PLUS_D3 = "r2scan_plus_d3"
        PBE = "pbe"
        PBE_PLUS_D3 = "pbe_plus_d3"
        PBE_U = "pbe_u"
        PBE_U_PLUS_D3 = "pbe_u_plus_d3"
        WB97XD = "wb97xd"
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorCalcMode.from_str)    @staticmethod
        def from_str(label: str) -> "EstimatorCalcMode":
            return EstimatorCalcMode(label.lower())
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorMethodType)class EstimatorMethodType(Enum):
        """
        Enum class which is used to determine `method_type` in `Estimator`.
    
        :ivar PFVM: Inference is performed in pfvm mode. Inference of dftd3 is performed in torch mode.
        :ivar PFVM_D3_PFVM: Inferences of pfp and dftd3 are performed in pfvm mode.
        :ivar MNCORE: Inference is performed on MN-Core.
        :ivar AUTO: The inference method is automatically selected based on the system size.
        """
    
        PFVM = "pfvm"
        MNCORE = "mncore"
        PFVM_D3_PFVM = "pfvm_d3_pfvm"
        AUTO = "auto"
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorMethodType.from_str)    @staticmethod
        def from_str(label: str) -> "EstimatorMethodType":
            return EstimatorMethodType(label.lower())
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorElementStatus)class EstimatorElementStatus(Enum):
        """
        Enum class which is used to determine `element_status` in `Estimator`.
    
        :ivar Expected: Fully supported elements.
        :ivar Experimental: Experimentally supported elements.
        """
    
        Expected = "EstimatorElementStatus.Expected"
        Experimental = "EstimatorElementStatus.Experimental"
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.unsupported_calc_mode_error)def unsupported_calc_mode_error(_calc_mode: Optional[Union[str, EstimatorCalcMode]]) -> None:
        _quote = "'"
        raise ValueError(
            f"{_calc_mode} is not a valid calculation mode. "
            f"Valid calculation modes are "
            f"{', '.join(_quote + str(mode.name) + _quote for mode in EstimatorCalcMode)}."
        )
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.unsupported_method_type_error)def unsupported_method_type_error(_method_type: Optional[Union[str, EstimatorMethodType]]) -> None:
        _quote = "'"
        raise ValueError(
            f"{_method_type} is not a valid inference method type. "
            f"Valid method types are "
            f"{', '.join(_quote + str(mode.name) + _quote for mode in EstimatorMethodType)}."
        )
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.get_calc_mode_enum)def get_calc_mode_enum(
        model_version: str,
        calc_mode: Optional[Union[EstimatorCalcMode, str]],
    ) -> Optional[EstimatorCalcMode]:
        _calc_mode = calc_mode.lower() if isinstance(calc_mode, str) else calc_mode
    
        if _calc_mode not in _calc_mode_table:
            unsupported_calc_mode_error(_calc_mode)
    
        # Check if _calc_mode is already of type EstimatorCalcMode
        if isinstance(_calc_mode, EstimatorCalcMode):
            return _calc_mode
    
        # Check if _calc_mode is a string and convert it
        if isinstance(_calc_mode, str):
            return EstimatorCalcMode.from_str(_calc_mode)
    
        # Handle the case where _calc_mode is None
        if _calc_mode is None:
            env_calc_mode = os.environ.get("MATLANTIS_PFP_CALC_MODE", None)
            if env_calc_mode is not None:
                if env_calc_mode.lower() in _calc_mode_table:
                    return EstimatorCalcMode.from_str(env_calc_mode)
                else:
                    unsupported_calc_mode_error(env_calc_mode)
            return _default_calc_mode_table_by_model_version.get(model_version)
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.get_method_type_enum)def get_method_type_enum(
        method_type: Optional[Union[EstimatorMethodType, str]],
    ) -> Optional[EstimatorMethodType]:
        _method_type = method_type.lower() if isinstance(method_type, str) else method_type
    
        if _method_type not in _method_type_table:
            unsupported_method_type_error(_method_type)
    
        # Check if _method_type is already of type EstimatorMethodType
        if isinstance(_method_type, EstimatorMethodType):
            return _method_type
    
        # Check if _method_type is a string and convert it
        if isinstance(_method_type, str):
            return EstimatorMethodType.from_str(_method_type)
    
        # Handle the case where _method_type is None
        if _method_type is None:
            env_method_type = os.environ.get("MATLANTIS_PFP_METHOD_TYPE") or os.environ.get(
                "PFP_INFERENCE_METHOD_TYPE"
            )
            if env_method_type is not None:
                if env_method_type.lower() in _method_type_table:
                    return EstimatorMethodType.from_str(env_method_type)
                else:
                    unsupported_method_type_error(env_method_type)
            return EstimatorMethodType.from_str(Config.method_type)
    
    
    
    
    _calc_mode_table: Mapping[Union[None, EstimatorCalcMode, str], Optional["CalcMode.ValueType"]] = {
        None: CALC_MODE_UNKNOWN,
        EstimatorCalcMode.CRYSTAL: CALC_MODE_CRYSTAL,
        EstimatorCalcMode.CRYSTAL.name: CALC_MODE_CRYSTAL,
        EstimatorCalcMode.CRYSTAL.value: CALC_MODE_CRYSTAL,
        EstimatorCalcMode.CRYSTAL_PLUS_D3: CALC_MODE_CRYSTAL_PLUS_D3,
        EstimatorCalcMode.CRYSTAL_PLUS_D3.name: CALC_MODE_CRYSTAL_PLUS_D3,
        EstimatorCalcMode.CRYSTAL_PLUS_D3.value: CALC_MODE_CRYSTAL_PLUS_D3,
        EstimatorCalcMode.CRYSTAL_U0: CALC_MODE_CRYSTAL_U0,
        EstimatorCalcMode.CRYSTAL_U0.name: CALC_MODE_CRYSTAL_U0,
        EstimatorCalcMode.CRYSTAL_U0.value: CALC_MODE_CRYSTAL_U0,
        EstimatorCalcMode.CRYSTAL_U0_PLUS_D3: CALC_MODE_CRYSTAL_U0_PLUS_D3,
        EstimatorCalcMode.CRYSTAL_U0_PLUS_D3.name: CALC_MODE_CRYSTAL_U0_PLUS_D3,
        EstimatorCalcMode.CRYSTAL_U0_PLUS_D3.value: CALC_MODE_CRYSTAL_U0_PLUS_D3,
        EstimatorCalcMode.MOLECULE: CALC_MODE_MOLECULE,
        EstimatorCalcMode.MOLECULE.name: CALC_MODE_MOLECULE,
        EstimatorCalcMode.MOLECULE.value: CALC_MODE_MOLECULE,
        EstimatorCalcMode.R2SCAN: CALC_MODE_R2SCAN,
        EstimatorCalcMode.R2SCAN.name: CALC_MODE_R2SCAN,
        EstimatorCalcMode.R2SCAN.value: CALC_MODE_R2SCAN,
        EstimatorCalcMode.R2SCAN_PLUS_D3: CALC_MODE_R2SCAN_PLUS_D3,
        EstimatorCalcMode.R2SCAN_PLUS_D3.name: CALC_MODE_R2SCAN_PLUS_D3,
        EstimatorCalcMode.R2SCAN_PLUS_D3.value: CALC_MODE_R2SCAN_PLUS_D3,
        EstimatorCalcMode.PBE: CALC_MODE_PBE,
        EstimatorCalcMode.PBE.name: CALC_MODE_PBE,
        EstimatorCalcMode.PBE.value: CALC_MODE_PBE,
        EstimatorCalcMode.PBE_PLUS_D3: CALC_MODE_PBE_PLUS_D3,
        EstimatorCalcMode.PBE_PLUS_D3.name: CALC_MODE_PBE_PLUS_D3,
        EstimatorCalcMode.PBE_PLUS_D3.value: CALC_MODE_PBE_PLUS_D3,
        EstimatorCalcMode.PBE_U: CALC_MODE_PBE_U,
        EstimatorCalcMode.PBE_U.name: CALC_MODE_PBE_U,
        EstimatorCalcMode.PBE_U.value: CALC_MODE_PBE_U,
        EstimatorCalcMode.PBE_U_PLUS_D3: CALC_MODE_PBE_U_PLUS_D3,
        EstimatorCalcMode.PBE_U_PLUS_D3.name: CALC_MODE_PBE_U_PLUS_D3,
        EstimatorCalcMode.PBE_U_PLUS_D3.value: CALC_MODE_PBE_U_PLUS_D3,
        EstimatorCalcMode.WB97XD: CALC_MODE_WB97XD,
        EstimatorCalcMode.WB97XD.name: CALC_MODE_WB97XD,
        EstimatorCalcMode.WB97XD.value: CALC_MODE_WB97XD,
    }
    
    _method_type_table: Mapping[
        Union[None, EstimatorMethodType, str], Optional["MethodType.ValueType"]
    ] = {
        None: None,
        EstimatorMethodType.PFVM: METHOD_TYPE_PFVM,
        EstimatorMethodType.PFVM.name: METHOD_TYPE_PFVM,
        EstimatorMethodType.PFVM.value: METHOD_TYPE_PFVM,
        EstimatorMethodType.MNCORE: METHOD_TYPE_MNCORE,
        EstimatorMethodType.MNCORE.name: METHOD_TYPE_MNCORE,
        EstimatorMethodType.MNCORE.value: METHOD_TYPE_MNCORE,
        EstimatorMethodType.PFVM_D3_PFVM: METHOD_TYPE_PFVM_D3_PFVM,
        EstimatorMethodType.PFVM_D3_PFVM.name: METHOD_TYPE_PFVM_D3_PFVM,
        EstimatorMethodType.PFVM_D3_PFVM.value: METHOD_TYPE_PFVM_D3_PFVM,
        EstimatorMethodType.AUTO: METHOD_TYPE_AUTO,
        EstimatorMethodType.AUTO.name: METHOD_TYPE_AUTO,
        EstimatorMethodType.AUTO.value: METHOD_TYPE_AUTO,
    }
    
    _default_calc_mode_table_by_model_version: Mapping[str, EstimatorCalcMode] = {
        "latest": EstimatorCalcMode.PBE,
        "v8.0.0": EstimatorCalcMode.PBE,
        "v7.0.0": EstimatorCalcMode.PBE,
        "v6.0.0": EstimatorCalcMode.PBE,
        "v5.0.0": EstimatorCalcMode.PBE,
        "v4.0.0": EstimatorCalcMode.PBE,
        "v3.0.0": EstimatorCalcMode.PBE_U,
        "v2.0.0": EstimatorCalcMode.PBE_U,
        "v0.0.0": EstimatorCalcMode.PBE_U,
    }
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorSystem)@dataclass(init=False)
    class EstimatorSystem:
        """
        Input argument of `Estimator.estimate` function.
    
        Parameters
        --------
        properties: List[str]
            What properties should be calculated. Current implemented properties are "energy", \
            "forces", "charges", and "virial". If gradient-based parameters (forces, virial) are not \
            requested, the estimator does not calculate them.
        coordinates: np.ndarray (dtype=np.float32 or np.float64)
            Coordination of atoms. The shape is `(n_atoms, 3)`.
        cell: np.ndarray (dtype=np.float32 or np.float64)
            Simulation cell of the system. The shape is `(3, 3)`
        species: Optional[np.ndarray] (dtype="<U3")
            Element list of atoms as string. Either species or atomic_numbers should be provided. \
            The shape is `(n_atoms, )`
        atomic_numbers: Optional[np.ndarray] (dtype=np.uint8)
            Element list of atoms as atomic number. Either species or atomic_numbers should be \
            provided. The shape is `(n_atoms, )`
        pbc: Optional[np.ndarray] (dtype=np.uint8)
            Whether periodic boundary conditions are applied (1) or not (0) for each axis. The shape is `(3, )`
        calc_mode: Optional[EstimatorCalcMode]
            Calculate results with given calculation mode.
        """
    
        properties: List[str]
        coordinates: np.ndarray
        cell: np.ndarray
        species: Optional[np.ndarray] = None
        atomic_numbers: Optional[np.ndarray] = None
        pbc: Optional[np.ndarray] = None
        calc_mode: Optional[EstimatorCalcMode] = None
    
        def __init__(
            self,
            properties: List[str],
            coordinates: np.ndarray,
            cell: np.ndarray,
            species: Optional[np.ndarray] = None,
            atomic_numbers: Optional[np.ndarray] = None,
            pbc: Optional[np.ndarray] = None,
            calc_mode: Optional[EstimatorCalcMode] = None,
            method_type: Optional[EstimatorMethodType] = None,
        ) -> None:
            self.properties = properties
            self.coordinates = coordinates
            self.cell = cell
            if isinstance(pbc, np.ndarray):
                pbc = pbc.astype(np.uint8)
            if isinstance(species, np.ndarray):
                species = species.astype("<U3")
            elif isinstance(atomic_numbers, np.ndarray):
                atomic_numbers = atomic_numbers.astype(np.uint8)
            self.pbc = pbc
            self.species = species
            self.atomic_numbers = atomic_numbers
            self.calc_mode = calc_mode
            self.method_type = method_type
    
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator)class Estimator:
        """
        Core class of PFP. Estimator calculates the energy and related values (e.g. atomic forces) \
        of a given atomic structure. A single `Estimator` instance can only be used by at most
        one `ASECalculator` instance.
    
        Acceptable PFP model versions are:
    
        - `v8.0.0`: latest version, the model with support for r2SCAN and its D3 correction.
    
        - `v7.0.0`: supporting 24 additional elements, for a total of 96 elements.
    
        - `v6.0.0`: improving the reproducibility of sparse environments
          (e.g. molecules in vacuum) and high-density structures.
    
        - `v5.0.0`: improving the reproducibility of physical properties
          for partially element-substituted structures, liquid, and organic crystals.
    
        - `v4.0.0`: improving the reproducibility of the energy order
          and the density of bulk structures of crystals, and intermolecular interactions.
    
        - `v3.0.0`: supporting 17 additional elements, for a total of 72 elements.
    
        - `v2.0.0`: model version with support for D3 correction and zero Hubbard U parameter.
    
        - `v0.0.0`: PFP paper version, does not support D3.
    
        - `latest`: alias for `v8.0.0`.
        """
    
        def __init__(
            self,
            max_retries: int = 10,
            model_version: Optional[str] = None,
            address: str = Config.address,
            priority: Optional[int] = None,
            calc_mode: Optional[Union[EstimatorCalcMode, str]] = None,
            method_type: Optional[Union[EstimatorMethodType, str]] = None,
        ) -> None:
            """
            Parameters
            --------
            max_retries: int, optional
                Maximum number of retries for requests. Default is 10.
            model_version: str, optional
                The model version. If not specified, the latest version will be used.
            address: str, optional
                The address of the PFP server. On Matlantis environment,
                users don't need to specify since this is automatically prepared.
            priority: int, optional
                Priority of the request. If not specified, Config.priority will be used.
            calc_mode: Optional[Union[EstimatorCalcMode, str]], optional
                Calculation mode. If not specified, the default calc_mode for the model_version will be used.
            method_type: Optional[Union[EstimatorMethodType, str]], optional
                Inference method type. If not specified, the default method_type from Config will be used.
            """
            self.channel = ChannelSingleton.get_channel(address=address)
            if not isinstance(max_retries, int):
                raise TypeError(
                    f'Expected type "int" for max_retries, got "{max_retries.__class__.__name__}"'
                )
            self.channel = grpc.intercept_channel(
                self.channel, RetryClientInterceptor(address, max_retries=max_retries)
            )
            # ErrorHandlerClientInterceptor must be placed at the last.
            self.channel = grpc.intercept_channel(self.channel, ErrorHandlerClientInterceptor(address))
            self.stub = pfp_pb2_grpc.EstimatorStub(self.channel)  # type: ignore
            # make it a Sequence so that nobody can modify the list accidentally
            self.implemented_properties: Sequence[str] = [
                "energy",
                "forces",
                "charges",
                "virial",
            ]
    
            if model_version is None:
                model_version = os.environ.get("MATLANTIS_PFP_MODEL_VERSION", latest_version)
            self._model_version: str = latest_version if model_version == "latest" else model_version
    
            self._calc_mode = get_calc_mode_enum(model_version=model_version, calc_mode=calc_mode)
            self._method_type = get_method_type_enum(method_type)
            self.message_isactive: Dict[MessageEnum, bool] = {m: True for m in MessageEnum}
    
            self._default_metadata: List[Tuple[str, str]] = [
                ("client-process-id", str(os.getpid())),
            ]
    
            self.model_api_stub = pfp_pb2_grpc.ModelAPIStub(self.channel)  # type: ignore
            metadata = self._default_metadata.copy()
            metadata.append((METADATA_METHOD_TYPE, str(self.method_type.value)))
            self.deployed_models: List[str] = self.model_api_stub.GetDeployedModels(
                pfp_pb2.ModelAPIRequest(),
                metadata=metadata,
            ).models
    
            if self._model_version not in self.deployed_models:
                raise ValueError(
                    f"Unexpected model version {self._model_version} was specified. "
                    f"Available model versions are {', '.join(self.deployed_models)}."
                )
    
            if priority is None:
                self.priority = Config.priority
            else:
                self.priority = priority
            self._application_context: Optional[str] = None
    
            self._attach_lock = threading.Lock()
            self._is_attached_to_calculator = False
            logger.debug(
                f"Created the estimator(model_version={self.model_version}, calc_mode={self.calc_mode}, method_type={self.method_type})"
            )
    
            self._path_table = {}  # type: ignore
            for _model_version in self.deployed_models:
                if _model_version == "latest":
                    continue
                self._path_table[_model_version] = {
                    CALC_MODE_MOLECULE: gaussian_shift_energies,
                    CALC_MODE_WB97XD: gaussian_shift_energies,
                }  # type: ignore
                _pt = self._path_table[_model_version]
                if version.parse(_model_version) >= version.parse("v9.0.0alpha3"):
                    _pt[CALC_MODE_CRYSTAL] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_CRYSTAL_U0] = vasp_wo_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE_U] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE] = vasp_wo_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_R2SCAN] = r2scan_wo_U_shift_energies_v1_9_3
                    _pt[CALC_MODE_R2SCAN_PLUS_D3] = _pt[CALC_MODE_R2SCAN]
                elif version.parse(_model_version) >= version.parse("v8.0.0"):
                    _pt[CALC_MODE_CRYSTAL] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_CRYSTAL_U0] = vasp_wo_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE_U] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE] = vasp_wo_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_R2SCAN] = r2scan_wo_U_shift_energies_v1_8_0
                    _pt[CALC_MODE_R2SCAN_PLUS_D3] = _pt[CALC_MODE_R2SCAN]
                elif version.parse(_model_version) >= version.parse("v7.0.0"):
                    _pt[CALC_MODE_CRYSTAL] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_CRYSTAL_U0] = vasp_wo_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE_U] = vasp_U_shift_energies_v1_7_1
                    _pt[CALC_MODE_PBE] = vasp_wo_U_shift_energies_v1_7_1
                elif version.parse(_model_version) >= version.parse("v4.0.0"):
                    _pt[CALC_MODE_CRYSTAL] = vasp_U_shift_energies_v1_4_0
                    _pt[CALC_MODE_CRYSTAL_U0] = vasp_wo_U_shift_energies_v1_4_0
                    _pt[CALC_MODE_PBE_U] = vasp_U_shift_energies_v1_4_0
                    _pt[CALC_MODE_PBE] = vasp_wo_U_shift_energies_v1_4_0
                else:
                    _pt[CALC_MODE_CRYSTAL] = vasp_U_shift_energies_v1_0_0
                    _pt[CALC_MODE_CRYSTAL_U0] = vasp_wo_U_shift_energies_v1_0_0
                    _pt[CALC_MODE_PBE_U] = vasp_U_shift_energies_v1_0_0
                    _pt[CALC_MODE_PBE] = vasp_wo_U_shift_energies_v1_0_0
                _pt[CALC_MODE_CRYSTAL_PLUS_D3] = _pt[CALC_MODE_CRYSTAL]
                _pt[CALC_MODE_CRYSTAL_U0_PLUS_D3] = _pt[CALC_MODE_CRYSTAL_U0]
                _pt[CALC_MODE_PBE_U_PLUS_D3] = _pt[CALC_MODE_CRYSTAL]
                _pt[CALC_MODE_PBE_PLUS_D3] = _pt[CALC_MODE_CRYSTAL_U0]
    
        def __repr__(self) -> str:
            return "<{} object; model_version={}, calc_mode={}, method_type={}>".format(
                type(self), self.model_version, self.calc_mode, self.method_type
            )
    
        def __enter__(self) -> "Estimator":
            return self
    
        def __exit__(
            self,
            ex_type: Optional[Type[Exception]],
            ex_value: Optional[Exception],
            trace: Optional[TracebackType],
        ) -> None:
            self.close()
    
        def close(self) -> None:
            """"""
            # NOTE: Do not close channel because of using shared channel
            # cf. https://github.pfidev.jp/Matlantis/pfp-api-client/pull/96
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.attach_to_calculator)    def attach_to_calculator(self, _: Calculator) -> None:
            with self._attach_lock:
                if self._is_attached_to_calculator:
                    raise MultiCalculatorUseDetected()
                # do not actually keep a reference to prevent circular references; we
                # could use a weakref if needed to keep a link to the owning Calculator
                self._is_attached_to_calculator = True
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.set_message_status)    def set_message_status(self, message: MessageEnum, message_enable: bool) -> None:
            """
            Toggle whether to show individual messages in `estimate()`.
    
            Parameters
            --------
            message: MessageEnum
                Message to toggle.
            message_enable: bool
                True/False correspond to enabling/disabling the message.
            """
            self.message_isactive[message] = message_enable
    
    
    
        @property
        def calc_mode(self) -> Optional[EstimatorCalcMode]:
            """The current calculation mode."""
            return self._calc_mode
    
        @calc_mode.setter
        def calc_mode(self, calc_mode: Union[EstimatorCalcMode, str]) -> None:
            """
            Set estimator calculation mode.
    
            Parameters
            --------
            calc_mode: Union[EstimatorCalcMode, str]
                The calculation mode to be set.
            """
    
            _calc_mode = calc_mode.lower() if isinstance(calc_mode, str) else calc_mode
    
            if _calc_mode not in _calc_mode_table:
                _quote = "'"
                raise ValueError(
                    f"{_calc_mode} is not a valid calculation mode. "
                    f"Valid calculation modes are "
                    f"{', '.join(_quote + str(mode.name) + _quote for mode in EstimatorCalcMode)}."
                )
    
            if isinstance(_calc_mode, str):
                _calc_mode = EstimatorCalcMode.from_str(_calc_mode)
    
            self._calc_mode = _calc_mode  # type: ignore
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.set_calc_mode)    def set_calc_mode(self, calc_mode: Union[EstimatorCalcMode, str]) -> None:
            """
            Set estimator calculation mode.
    
            Parameters
            --------
            calc_mode: Union[EstimatorCalcMode, str]
                The calculation mode to be set.
            """
            self.calc_mode = calc_mode  # type: ignore
    
    
    
        @property
        def method_type(self) -> EstimatorMethodType:
            return self._method_type  # type: ignore
    
        @method_type.setter
        def method_type(self, method_type: Union[EstimatorMethodType, str]) -> None:
            """
            Set estimator method type.
    
            Parameters
            --------
            method_type: Union[EstimatorMethodType, str]
                The inference method type to be set.
            """
    
            _method_type: Union[EstimatorMethodType, str]
            _method_type = method_type.lower() if isinstance(method_type, str) else method_type
    
            if _method_type not in _method_type_table:
                _quote = "'"
                raise ValueError(
                    f"{_method_type} is not a valid inference method type. "
                    f"Valid method types are "
                    f"{', '.join(_quote + str(mode.name) + _quote for mode in EstimatorMethodType)}."
                )
            if isinstance(_method_type, str):
                _method_type = EstimatorMethodType.from_str(_method_type)
            self._method_type = _method_type
            metadata = self._default_metadata.copy()
            metadata.append((METADATA_METHOD_TYPE, str(self.method_type.value)))
            self.deployed_models = self.model_api_stub.GetDeployedModels(
                pfp_pb2.ModelAPIRequest(),
                metadata=metadata,
            ).models
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.set_method_type)    def set_method_type(self, method_type: Union[EstimatorMethodType, str]) -> None:
            """
            Set estimator method type.
    
            Parameters
            --------
            method_type: Union[EstimatorMethodType, str]
                The inference method type to be set.
            """
            self.method_type = method_type  # type: ignore
    
    
    
        @property
        def available_models(self) -> Sequence[str]:
            """Return the list of model versions available on the server."""
            return self.deployed_models
    
        @property
        def model_version(self) -> str:
            """The current model version."""
            return self._model_version
    
        @model_version.setter
        def model_version(self, model_version: str) -> None:
            """
            Set model version.
    
            Parameters
            --------
            model_version: str
                The model version to be set.
            """
            if (
                model_version == "latest" and self._model_version == latest_version
            ) or self._model_version == model_version:
                return
    
            if model_version not in self.deployed_models:
                raise ValueError(
                    f"Unexpected model version was specified. Available model versions are {', '.join(self.deployed_models)}."
                )
            self._model_version = latest_version if model_version == "latest" else model_version
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.set_model_version)    def set_model_version(self, model_version: str) -> None:
            """
            Set model version.
    
            Parameters
            --------
            model_version: str
                The model version to be set.
            """
            self.model_version = model_version
    
    
    
        @property
        def priority(self) -> int:
            """The current request priority."""
            return self._priority
    
        @priority.setter
        def priority(self, priority: int) -> None:
            """
            Set priority of the request.
    
            Parameters
            --------
            priority: int
                Priority of the request. Must be within 1-100.
            """
            if not isinstance(priority, int):
                raise ValueError(f"priority must be int, but {type(priority)}")
            if priority <= 0 or 100 < priority:
                raise ValueError(f"priority must be within 1-100, but {priority}.")
            self._priority = priority
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.set_priority)    def set_priority(self, priority: int) -> None:
            """
            Set priority of the request.
    
            Parameters
            --------
            priority: int
                Priority of the request.
            """
            self.priority = priority
    
    
    
        @property
        def application_context(self) -> Optional[str]:
            return self._application_context
    
        @application_context.setter
        def application_context(self, application_context: Optional[str]) -> None:
            self._application_context = application_context
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.supported_elements)    def supported_elements(
            self,
            model_version: Optional[str] = None,
            status: Union[
                EstimatorElementStatus, Tuple[EstimatorElementStatus, ...]
            ] = EstimatorElementStatus.Expected,
            calc_mode: Optional[Union[EstimatorCalcMode, str]] = None,
        ) -> List[str]:
            """
            Return list of implemented elements in specified model version.
    
            Parameters
            --------
            model_version: str, optional
                Model version to be specified.
            status: EstimatorElementStatus or tuple of EstimatorElementStatus, optional
                Target support status to return. Multiple keywords can be specified.
            calc_mode: EstimatorCalcMode, or str, optional
                Calculation mode. If not specified, calc_mode set in the estimator will be used.
    
            Returns
            --------
            elements: List[str]
                List which contains implemented elements.
            """
    
            _status_table: Dict[EstimatorElementStatus, "TargetElementStatus.ValueType"] = {
                EstimatorElementStatus.Expected: ELEMENT_STATUS_EXPECTED,
                EstimatorElementStatus.Experimental: ELEMENT_STATUS_EXPERIMENTAL,
            }
    
            if model_version == "latest":
                model_version = latest_version
    
            elif model_version is None:
                model_version = self.model_version
    
            if model_version not in self.deployed_models:
                raise ValueError(
                    f"Unexpected model version was specified.Available model versions are {', '.join(self.deployed_models)}."
                )
    
            _target = set((status,)) if isinstance(status, EstimatorElementStatus) else set(status)
            _target_element_status: Set["TargetElementStatus.ValueType"] = set()
    
            for t in _target:
                if t not in _status_table:
                    raise ValueError(
                        f"Unexpected status was specified. Available status are {', '.join([s.value for s in _status_table.keys()])}."
                    )
                _target_element_status.add(_status_table[t])
    
            _calc_mode: Optional[Union[EstimatorCalcMode, str]]
    
            if calc_mode is not None:
                _calc_mode = calc_mode.lower() if isinstance(calc_mode, str) else calc_mode
            else:
                _calc_mode = self._calc_mode
    
            if _calc_mode in _calc_mode_table:
                request = pfp_pb2.GetSupportedElementsRequest(
                    model_version=model_version,
                    target_element_status=_target_element_status,
                    calc_mode=_calc_mode_table[_calc_mode],  # type: ignore[arg-type]
                )
                elements: List[str] = list(self.model_api_stub.GetSupportedElements(request).elements)
                return elements
            else:
                raise ValueError(
                    f"{_calc_mode} is not a valid calculation mode. "
                    f"Valid calculation modes are {', '.join(mode.name for mode in EstimatorCalcMode)}."
                )
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.get_shift_energy_table)    def get_shift_energy_table(
            self,
            calc_mode: Union[None, EstimatorCalcMode, str] = None,
            model_version: Optional[str] = None,
        ) -> Dict[int, float]:
            """
            Return the table of shift energies of input atoms.
    
            Parameters
            --------
            calc_mode: str, EstimatorCalcMode, or None, optional
                Calculation mode. If not specified, calc_mode set in the estimator will be used.
            model_version: str, optional
                The model version. If not specified, model_version set in the estimator will be used.
    
            Returns
            --------
            shift_energies_table: Dict[int, float]
                Table of shift energies.
            """
    
            if calc_mode is None:
                calc_mode = self._calc_mode
    
            if calc_mode is None:
                raise ValueError(
                    "Please set calc_mode in the estimator or specify calc_mode explicitly."
                )
    
            _calc_mode = calc_mode.lower() if isinstance(calc_mode, str) else calc_mode
    
            if _calc_mode not in _calc_mode_table:
                _quote = "'"
                raise ValueError(
                    f"{_calc_mode} is not a valid calculation mode. "
                    f"Valid calculation modes are "
                    f"{', '.join(_quote + str(mode.name) + _quote for mode in EstimatorCalcMode)}."
                )
    
            if model_version == "latest":
                model_version = latest_version
    
            elif model_version is None:
                model_version = self.model_version
    
            if model_version not in self.deployed_models:
                raise ValueError(
                    f"Unexpected model version was specified.Available model versions are {', '.join(self.deployed_models)}."
                )
            if _calc_mode in [
                EstimatorCalcMode.R2SCAN,
                EstimatorCalcMode.R2SCAN.name,
                EstimatorCalcMode.R2SCAN.value,
            ] and version.parse(model_version) < version.parse("v8.0.0"):
                raise ValueError("The r2scan calc_mode is supported by v8.0.0 and newer models.")
    
            _target = self._path_table[model_version][_calc_mode_table[_calc_mode]]  # type: ignore
            shift_energies_table: Dict[int, float] = _target.shift_energies  # type: ignore
            for i in range(1, 128):  # Unsupported elements are filled with 0
                if i not in shift_energies_table:
                    shift_energies_table[i] = 0.0
    
            return shift_energies_table
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.estimate)    def estimate(
            self, args: EstimatorSystem
        ) -> Dict[str, Union[List[MessageEnum], np.ndarray, Dict[str, int]]]:
            """
            This function receives atomic system information and returns energy, forces and related \
            values.
    
            Parameters
            --------
            args: EstimatorSystem
                For details, see the document for EstimatorSystem.
    
            Returns
            --------
            results: Dict[str, any]
                Dictionary which contains various calculated properties
            """
    
            if not args.properties:
                args.properties = list(self.implemented_properties)
    
            if (
                self.model_version == "latest"
                and max(map(version.parse, self.available_models)) >= version.parse("v2.0.0")
            ) or version.parse(self.model_version) >= version.parse("v2.0.0"):
                self.coordinates: np.ndarray = args.coordinates.astype(np.float64)
                self.cell: np.ndarray = args.cell.astype(np.float64)
            else:
                self.coordinates = args.coordinates.astype(np.float32)
                self.cell = args.cell.astype(np.float32)
    
            coord_ndarray_pb = to_ndarray_pb(self.coordinates)  # np.float32 or np.float64
            cell_ndarray_pb = to_ndarray_pb(self.cell)  # np.float32 or np.float64
            atomic_numbers_ndarray_pb: Optional[NdArray]
            if isinstance(args.atomic_numbers, np.ndarray):
                # np.uint8
                atomic_numbers_ndarray_pb = to_ndarray_pb(args.atomic_numbers)
                species_ndarray_pb = None
            elif isinstance(args.species, np.ndarray):
                atomic_numbers_ndarray_pb = None
                species_ndarray_pb = to_ndarray_pb(args.species)  # np.uint8
            else:
                raise ValueError("Must specify either atomic_numbers or species")
            # np.uint8
            pbc_ndarray_pb = to_ndarray_pb(args.pbc)  # type: ignore
            _calc_mode: Optional[EstimatorCalcMode]
            if args.calc_mode is not None:
                _calc_mode = args.calc_mode
            else:
                _calc_mode = self._calc_mode
            _method_type = self._method_type
            if _calc_mode in _calc_mode_table:
                request = pfp_pb2.EstimateRequest(
                    model_version=self.model_version,
                    properties=args.properties,
                    atomic_numbers=atomic_numbers_ndarray_pb,
                    species=species_ndarray_pb,
                    cell=cell_ndarray_pb,
                    pbc=pbc_ndarray_pb,
                    coordinates=coord_ndarray_pb,
                    calc_mode=_calc_mode_table[_calc_mode],  # type: ignore
                    method_type=_method_type_table[_method_type],  # type: ignore
                )
                metadata = self._default_metadata.copy()
                metadata.append((METADATA_N_ATOMS, str(len(args.coordinates))))
                metadata.append((METADATA_MODEL_VERSION, str(self.model_version)))
                metadata.append((METADATA_CALC_MODE, str(_calc_mode_table[_calc_mode])))
                metadata.append((METADATA_METHOD_TYPE, str(self.method_type.value)))
                metadata.append((METADATA_CLIENT_SIDE_PRIORITY, str(self.priority)))
                metadata.append((METADATA_CLIENT_VERSION, __version__))
                metadata.append(
                    (
                        METADATA_REQUEST_ID,
                        f"{Config.identifier}-{int(time.time() * 1000)}-{os.getpid()}-{threading.get_ident()}",
                    )
                )
                if (execution_context := os.environ.get("MATLANTIS_EXECUTION_CONTEXT")) is not None:
                    metadata.append((METADATA_EXECUTION_CONTEXT, execution_context))
                application_context = (
                    self.application_context
                    if self.application_context is not None
                    else os.environ.get(
                        "MATLANTIS_APPLICATION_CONTEXT",
                        "pfp-no-context",
                    )
                )
                metadata.append((METADATA_APPLICATION_CONTEXT, application_context.lower()))
                response = self.stub.Estimate(request, metadata=metadata)
            else:
                raise ValueError(
                    f"{_calc_mode} is not a valid calculation mode. "
                    f"Valid calculation modes are {', '.join(mode.name for mode in EstimatorCalcMode)}."
                )
    
            results = {}
            results["energy"] = response.energy
            if "forces" in args.properties:
                results["forces"] = from_ndarray_pb(response.forces)
            if "charges" in args.properties:
                # Reshape from (n_atoms, 1) to (n_atoms,)
                results["charges"] = from_ndarray_pb(response.charges).ravel()
            if "virial" in args.properties:
                results["virial"] = from_ndarray_pb(response.virial)
            results["messages"] = from_messages_pb(response.messages)
            results["calc_stats"] = dict(response.calc_stats)
            return results
    
    
    
    
    
    [[docs]](../../../pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator.get_descriptors)    def get_descriptors(
            self,
            args: EstimatorSystem,
            indices: Optional[np.ndarray] = None,
        ) -> Dict[str, Union[np.ndarray, List[MessageEnum], Dict[str, int]]]:
            """
            This function receives atomic system information and returns descriptors.
    
            Parameters
            --------
            args: EstimatorSystem
                For details, see the document for EstimatorSystem.
            indices: Optional[np.ndarray]
                A 1-dimensional array of indices (between 0 and n_atoms-1) specifying the
                atomic indices for descriptors to be obtained.
    
            Returns
            --------
            results: Dict[str, any]
                Dictionary which contains scalar descriptors and other calculation information.
                The `descriptors` key contains scalar descriptors with shape of [n_atoms, 256].
            """
            coordinates: np.ndarray = args.coordinates.astype(np.float64)
            cell: np.ndarray = args.cell.astype(np.float64)
    
            coord_ndarray_pb = to_ndarray_pb(coordinates)  # np.float64 only
            cell_ndarray_pb = to_ndarray_pb(cell)  # np.float64 only
            atomic_numbers_ndarray_pb: Optional[NdArray]
            if isinstance(args.atomic_numbers, np.ndarray):
                # np.uint8
                atomic_numbers_ndarray_pb = to_ndarray_pb(args.atomic_numbers)
                n_atoms = args.atomic_numbers.shape[0]
                species_ndarray_pb = None
            elif isinstance(args.species, np.ndarray):
                atomic_numbers_ndarray_pb = None
                species_ndarray_pb = to_ndarray_pb(args.species)  # np.uint8
                n_atoms = args.species.shape[0]
            else:
                raise ValueError("Must specify either atomic_numbers or species")
            # np.uint8
            pbc_ndarray_pb = to_ndarray_pb(args.pbc)  # type: ignore
            _calc_mode: Optional[EstimatorCalcMode]
            if args.calc_mode is not None:
                _calc_mode = args.calc_mode
            else:
                _calc_mode = self._calc_mode
            _method_type = self._method_type
            if _calc_mode in _calc_mode_table:
                if indices is not None:
                    if len(indices) == 0:
                        indices = None
                    else:
                        indices = indices.astype(np.uint32)
                        # there is also server-side validation for this
                        if not np.all((indices >= 0) & (indices < n_atoms)):
                            raise ValueError(
                                f"Indices for descriptors must be between 0 and {n_atoms - 1} (= n_atoms - 1)."
                            )
                        n_indices = len(indices)
                        indices = to_ndarray_pb(indices)  # type: ignore
    
                request = pfp_pb2.GetDescriptorsRequestV2(
                    model_version=self.model_version,
                    properties=args.properties,
                    atomic_numbers=atomic_numbers_ndarray_pb,
                    species=species_ndarray_pb,
                    cell=cell_ndarray_pb,
                    pbc=pbc_ndarray_pb,
                    coordinates=coord_ndarray_pb,
                    calc_mode=_calc_mode_table[_calc_mode],  # type: ignore
                    method_type=_method_type_table[_method_type],  # type: ignore
                    indices=indices,  # type: ignore
                )
                metadata = self._default_metadata.copy()
                metadata.append((METADATA_N_ATOMS, str(len(args.coordinates))))
                metadata.append((METADATA_MODEL_VERSION, str(self.model_version)))
                metadata.append((METADATA_CALC_MODE, str(_calc_mode_table[_calc_mode])))
                metadata.append((METADATA_METHOD_TYPE, str(self.method_type.value)))
                metadata.append((METADATA_CLIENT_SIDE_PRIORITY, str(self.priority)))
                metadata.append((METADATA_CLIENT_VERSION, __version__))
                metadata.append(
                    (
                        METADATA_REQUEST_ID,
                        f"{Config.identifier}-{int(time.time() * 1000)}-{os.getpid()}-{threading.get_ident()}",
                    )
                )
                if indices is not None:
                    metadata.append((METADATA_N_INDICES, str(n_indices)))
                if (execution_context := os.environ.get("MATLANTIS_EXECUTION_CONTEXT")) is not None:
                    metadata.append((METADATA_EXECUTION_CONTEXT, execution_context))
                application_context = (
                    self.application_context
                    if self.application_context is not None
                    else os.environ.get(
                        "MATLANTIS_APPLICATION_CONTEXT",
                        "pfp-no-context",
                    )
                )
                metadata.append((METADATA_APPLICATION_CONTEXT, application_context.lower()))
                response = self.stub.GetDescriptorsV2(request, metadata=metadata)
            else:
                raise ValueError(
                    f"{_calc_mode} is not a valid calculation mode. "
                    f"Valid calculation modes are {', '.join(mode.name for mode in EstimatorCalcMode)}."
                )
    
            results: Dict[str, Union[np.ndarray, List[MessageEnum], Dict[str, int]]] = {
                "descriptors": from_ndarray_pb(response.descriptors),
                "messages": from_messages_pb(response.messages),
                "calc_stats": cast(Dict[str, int], dict(response.calc_stats)),
            }
            return results
    
    
    
