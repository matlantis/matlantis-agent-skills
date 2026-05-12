# Estimator#

_class _pfp_api_client.pfp.estimator.Estimator(_max_retries : int = 10_, _model_version : Optional[str] = None_, _address : str = 'localhost:5000'_, _priority : Optional[int] = None_, _calc_mode : Optional[Union[str, EstimatorCalcMode]] = None_, _method_type : Optional[Union[str, EstimatorMethodType]] = None_)[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator)#
    

Bases: `object`

Core class of PFP. Estimator calculates the energy and related values (e.g. atomic forces) of a given atomic structure. A single Estimator instance can only be used by at most one ASECalculator instance.

Acceptable PFP model versions are:

  * v8.0.0: latest version, the model with support for r2SCAN and its D3 correction.

  * v7.0.0: supporting 24 additional elements, for a total of 96 elements.

  * v6.0.0: improving the reproducibility of sparse environments (e.g. molecules in vacuum) and high-density structures.

  * v5.0.0: improving the reproducibility of physical properties for partially element-substituted structures, liquid, and organic crystals.

  * v4.0.0: improving the reproducibility of the energy order and the density of bulk structures of crystals, and intermolecular interactions.

  * v3.0.0: supporting 17 additional elements, for a total of 72 elements.

  * v2.0.0: model version with support for D3 correction and zero Hubbard U parameter.

  * v0.0.0: PFP paper version, does not support D3.

  * latest: alias for v8.0.0.




_property _application_context _: Optional[str]_#
    

attach_to_calculator(__ : Calculator_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.attach_to_calculator)#
    

_property _available_models _: Sequence[str]_#
    

Return the list of model versions available on the server.

_property _calc_mode _: Optional[EstimatorCalcMode]_#
    

The current calculation mode.

estimate(_args : EstimatorSystem_) → Dict[str, Union[List[[MessageEnum](pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum "pfp_api_client.pfp.utils.messages.MessageEnum")], ndarray, Dict[str, int]]][[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.estimate)#
    

This function receives atomic system information and returns energy, forces and related values.

Parameters
    

**args** (_EstimatorSystem_) – For details, see the document for EstimatorSystem.

Returns
    

**results** – Dictionary which contains various calculated properties

Return type
    

Dict[str, any]

get_descriptors(_args : EstimatorSystem_, _indices : Optional[ndarray] = None_) → Dict[str, Union[List[[MessageEnum](pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum "pfp_api_client.pfp.utils.messages.MessageEnum")], ndarray, Dict[str, int]]][[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.get_descriptors)#
    

This function receives atomic system information and returns descriptors.

Parameters
    

  * **args** (_EstimatorSystem_) – For details, see the document for EstimatorSystem.

  * **indices** (_Optional_ _[__np.ndarray_ _]_) – A 1-dimensional array of indices (between 0 and n_atoms-1) specifying the atomic indices for descriptors to be obtained.



Returns
    

**results** – Dictionary which contains scalar descriptors and other calculation information. The descriptors key contains scalar descriptors with shape of [n_atoms, 256].

Return type
    

Dict[str, any]

get_shift_energy_table(_calc_mode : Union[None, EstimatorCalcMode, str] = None_, _model_version : Optional[str] = None_) → Dict[int, float][[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.get_shift_energy_table)#
    

Return the table of shift energies of input atoms.

Parameters
    

  * **calc_mode** (_str_ _,__EstimatorCalcMode_ _, or_ _None_ _,__optional_) – Calculation mode. If not specified, calc_mode set in the estimator will be used.

  * **model_version** (_str_ _,__optional_) – The model version. If not specified, model_version set in the estimator will be used.



Returns
    

**shift_energies_table** – Table of shift energies.

Return type
    

Dict[int, float]

_property _method_type _: EstimatorMethodType_#
    

_property _model_version _: str_#
    

The current model version.

_property _priority _: int_#
    

The current request priority.

set_calc_mode(_calc_mode : Union[EstimatorCalcMode, str]_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.set_calc_mode)#
    

Set estimator calculation mode.

Parameters
    

**calc_mode** (_Union_ _[__EstimatorCalcMode_ _,__str_ _]_) – The calculation mode to be set.

set_message_status(_message : [MessageEnum](pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum "pfp_api_client.pfp.utils.messages.MessageEnum")_, _message_enable : bool_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.set_message_status)#
    

Toggle whether to show individual messages in estimate().

Parameters
    

  * **message** ([_MessageEnum_](pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum "pfp_api_client.pfp.utils.messages.MessageEnum")) – Message to toggle.

  * **message_enable** (_bool_) – True/False correspond to enabling/disabling the message.




set_method_type(_method_type : Union[EstimatorMethodType, str]_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.set_method_type)#
    

Set estimator method type.

Parameters
    

**method_type** (_Union_ _[__EstimatorMethodType_ _,__str_ _]_) – The inference method type to be set.

set_model_version(_model_version : str_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.set_model_version)#
    

Set model version.

Parameters
    

**model_version** (_str_) – The model version to be set.

set_priority(_priority : int_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.set_priority)#
    

Set priority of the request.

Parameters
    

**priority** (_int_) – Priority of the request.

supported_elements(_model_version : Optional[str] = None_, _status : Union[EstimatorElementStatus, Tuple[EstimatorElementStatus, ...]] = EstimatorElementStatus.Expected_, _calc_mode : Optional[Union[str, EstimatorCalcMode]] = None_) → List[str][[source]](_modules/pfp_api_client/pfp/estimator.html#Estimator.supported_elements)#
    

Return list of implemented elements in specified model version.

Parameters
    

  * **model_version** (_str_ _,__optional_) – Model version to be specified.

  * **status** (_EstimatorElementStatus_ _or_ _tuple of EstimatorElementStatus_ _,__optional_) – Target support status to return. Multiple keywords can be specified.

  * **calc_mode** (_EstimatorCalcMode_ _, or_ _str_ _,__optional_) – Calculation mode. If not specified, calc_mode set in the estimator will be used.



Returns
    

**elements** – List which contains implemented elements.

Return type
    

List[str]

_class _pfp_api_client.pfp.estimator.EstimatorCalcMode(_value_ , _names =None_, _*_ , _module =None_, _qualname =None_, _type =None_, _start =1_, _boundary =None_)[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorCalcMode)#
    

Bases: `Enum`

Enum class which is used to determine calc_mode in Estimator.

Variables
    

  * **R2SCAN** – This mode corresponds to DFT calculations with r2SCAN functional and plane wave basis sets.

  * **R2SCAN_PLUS_D3** – R2SCAN mode with D3 correction enabled. This mode corresponds to DFT calculations with r2SCAN functional and plane wave basis sets, plus a D3 dispersion correction.

  * **PBE** – This mode corresponds to DFT calculations with PBE functional and plane wave basis sets.

  * **PBE_PLUS_D3** – PBE mode with D3 correction enabled. This mode corresponds to DFT calculations with plane wave basis sets, plus a D3 dispersion correction.

  * **PBE_U** – PBE mode with Hubbard U parameter. This mode corresponds to DFT calculations with plane wave basis sets.

  * **PBE_U_PLUS_D3** – PBE_U mode with D3 correction enabled. This mode corresponds to DFT calculations with plane wave basis sets, plus a D3 dispersion correction.

  * **WB97XD** – This mode corresponds to DFT calculations with local basis sets.




CRYSTAL _ = 'crystal'_#
    

CRYSTAL_PLUS_D3 _ = 'crystal_plus_d3'_#
    

CRYSTAL_U0 _ = 'crystal_u0'_#
    

CRYSTAL_U0_PLUS_D3 _ = 'crystal_u0_plus_d3'_#
    

MOLECULE _ = 'molecule'_#
    

PBE _ = 'pbe'_#
    

PBE_PLUS_D3 _ = 'pbe_plus_d3'_#
    

PBE_U _ = 'pbe_u'_#
    

PBE_U_PLUS_D3 _ = 'pbe_u_plus_d3'_#
    

R2SCAN _ = 'r2scan'_#
    

R2SCAN_PLUS_D3 _ = 'r2scan_plus_d3'_#
    

WB97XD _ = 'wb97xd'_#
    

_static _from_str(_label : str_) → EstimatorCalcMode[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorCalcMode.from_str)#
    

_class _pfp_api_client.pfp.estimator.EstimatorElementStatus(_value_ , _names =None_, _*_ , _module =None_, _qualname =None_, _type =None_, _start =1_, _boundary =None_)[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorElementStatus)#
    

Bases: `Enum`

Enum class which is used to determine element_status in Estimator.

Variables
    

  * **Expected** – Fully supported elements.

  * **Experimental** – Experimentally supported elements.




Expected _ = 'EstimatorElementStatus.Expected'_#
    

Experimental _ = 'EstimatorElementStatus.Experimental'_#
    

_class _pfp_api_client.pfp.estimator.EstimatorMethodType(_value_ , _names =None_, _*_ , _module =None_, _qualname =None_, _type =None_, _start =1_, _boundary =None_)[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorMethodType)#
    

Bases: `Enum`

Enum class which is used to determine method_type in Estimator.

Variables
    

  * **PFVM** – Inference is performed in pfvm mode. Inference of dftd3 is performed in torch mode.

  * **PFVM_D3_PFVM** – Inferences of pfp and dftd3 are performed in pfvm mode.

  * **MNCORE** – Inference is performed on MN-Core.

  * **AUTO** – The inference method is automatically selected based on the system size.




AUTO _ = 'auto'_#
    

MNCORE _ = 'mncore'_#
    

PFVM _ = 'pfvm'_#
    

PFVM_D3_PFVM _ = 'pfvm_d3_pfvm'_#
    

_static _from_str(_label : str_) → EstimatorMethodType[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorMethodType.from_str)#
    

_class _pfp_api_client.pfp.estimator.EstimatorSystem(_properties : List[str]_, _coordinates : ndarray_, _cell : ndarray_, _species : Optional[ndarray] = None_, _atomic_numbers : Optional[ndarray] = None_, _pbc : Optional[ndarray] = None_, _calc_mode : Optional[EstimatorCalcMode] = None_, _method_type : Optional[EstimatorMethodType] = None_)[[source]](_modules/pfp_api_client/pfp/estimator.html#EstimatorSystem)#
    

Bases: `object`

Input argument of Estimator.estimate function.

Parameters
    

  * **properties** (_List_ _[__str_ _]_) – What properties should be calculated. Current implemented properties are “energy”, “forces”, “charges”, and “virial”. If gradient-based parameters (forces, virial) are not requested, the estimator does not calculate them.

  * **coordinates** (_np.ndarray_ _(__dtype=np.float32_ _or_ _np.float64_ _)_) – Coordination of atoms. The shape is (n_atoms, 3).

  * **cell** (_np.ndarray_ _(__dtype=np.float32_ _or_ _np.float64_ _)_) – Simulation cell of the system. The shape is (3, 3)

  * **species** (_Optional_ _[__np.ndarray_ _]__(__dtype=" <U3"__)_) – Element list of atoms as string. Either species or atomic_numbers should be provided. The shape is (n_atoms, )

  * **atomic_numbers** (_Optional_ _[__np.ndarray_ _]__(__dtype=np.uint8_ _)_) – Element list of atoms as atomic number. Either species or atomic_numbers should be provided. The shape is (n_atoms, )

  * **pbc** (_Optional_ _[__np.ndarray_ _]__(__dtype=np.uint8_ _)_) – Whether periodic boundary conditions are applied (1) or not (0) for each axis. The shape is (3, )

  * **calc_mode** (_Optional_ _[__EstimatorCalcMode_ _]_) – Calculate results with given calculation mode.




atomic_numbers _: Optional[ndarray]__ = None_#
    

calc_mode _: Optional[EstimatorCalcMode]__ = None_#
    

cell _: ndarray_#
    

coordinates _: ndarray_#
    

pbc _: Optional[ndarray]__ = None_#
    

properties _: List[str]_#
    

species _: Optional[ndarray]__ = None_#
    

pfp_api_client.pfp.estimator.get_calc_mode_enum(_model_version : str_, _calc_mode : Optional[Union[str, EstimatorCalcMode]]_) → Optional[EstimatorCalcMode][[source]](_modules/pfp_api_client/pfp/estimator.html#get_calc_mode_enum)#
    

pfp_api_client.pfp.estimator.get_method_type_enum(_method_type : Optional[Union[str, EstimatorMethodType]]_) → Optional[EstimatorMethodType][[source]](_modules/pfp_api_client/pfp/estimator.html#get_method_type_enum)#
    

pfp_api_client.pfp.estimator.unsupported_calc_mode_error(__calc_mode : Optional[Union[str, EstimatorCalcMode]]_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#unsupported_calc_mode_error)#
    

pfp_api_client.pfp.estimator.unsupported_method_type_error(__method_type : Optional[Union[str, EstimatorMethodType]]_) → None[[source]](_modules/pfp_api_client/pfp/estimator.html#unsupported_method_type_error)#
    
