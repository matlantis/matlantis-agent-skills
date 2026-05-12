# Source code for matlantis_features.features.base
    
    
    import contextlib
    import functools
    import json
    import os
    import pathlib
    import tempfile
    import warnings
    from dataclasses import dataclass, fields
    from typing import Any, Callable, Dict, Iterator, List, Optional, Set, Type, TypeVar, Union
    
    from ase import Atoms
    from ga4mp import GtagMP
    
    # TODO minimize the dependency to pfp-api-client. cf. https://github.com/pfn-Matlantis/matlantis-features/pull/785
    from pfp_api_client.pfp.config import Config
    from pfp_api_client.pfp.estimator import _default_calc_mode_table_by_model_version
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.utils import Context, FeatureCost
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    from matlantis_features.utils.save_util import (
        NpEncoder,
        from_dict,
        get_version_info,
        to_serialized_dict,
    )
    
    T = TypeVar("T")
    Self = TypeVar("Self", bound="FeatureBase")
    
    
    def _is_atoms(var: Any) -> bool:
        return isinstance(var, Atoms) or isinstance(var, MatlantisAtoms)
    
    
    def _is_list_of_atoms(var: Any) -> bool:
        if isinstance(var, list):
            return all([_is_atoms(i) for i in var])
        return False
    
    
    def _is_list_of_list_of_atoms(var: Any) -> bool:
        if isinstance(var, list):
            return all([_is_list_of_atoms(i) for i in var])
        return False
    
    
    def _copy_atoms(atoms: Union[Atoms, MatlantisAtoms]) -> Union[Atoms, MatlantisAtoms]:
        if isinstance(atoms, Atoms):
            atoms_copy = atoms.copy()
            atoms_copy.calc = None
            return atoms_copy
        elif isinstance(atoms, MatlantisAtoms):
            ase_atoms_copy = atoms.ase_atoms.copy()
            ase_atoms_copy.calc = None
            atoms_copy = MatlantisAtoms(atoms=ase_atoms_copy)
            return atoms_copy
        raise TypeError("Input must be of type ASEAtoms or MatlantisAtoms.")
    
    
    def _copy_list_of_atoms(list_atoms: List[Union[Atoms, MatlantisAtoms]]) -> List[Union[Atoms, MatlantisAtoms]]:
        return [_copy_atoms(atoms) for atoms in list_atoms]
    
    
    def _copy_list_of_list_of_atoms(list_atoms: List[List[Union[Atoms, MatlantisAtoms]]]) -> List[List[Union[Atoms, MatlantisAtoms]]]:
        return [_copy_list_of_atoms(sublist) for sublist in list_atoms]
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.copy_to_child_conversion.html#matlantis_features.features.base.copy_to_child_conversion)def copy_to_child_conversion(var: Any) -> Any:
        """Convert the input variables of the __call__ function of the feature. If the \
        variable is ASEAtoms, the ASEAtoms will be copied and the attached calculator will be \
        removed.
    
        Args:
            var (Any): The input variables.
        Returns:
            Any : The variables after conversion.
        """
        if _is_atoms(var):
            return _copy_atoms(var)
        # copy List[Union[Atoms, MatlantisAtoms]] (e.g. NEBFeature)
        elif _is_list_of_atoms(var):
            return _copy_list_of_atoms(var)
        # copy List[List[Union[Atoms, MatlantisAtoms]]] (e.g. ReactionStringFeature)
        elif _is_list_of_list_of_atoms(var):
            return _copy_list_of_list_of_atoms(var)
        return var
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseCaller.html#matlantis_features.features.base.FeatureBaseCaller)class FeatureBaseCaller(type):
        """The metaclass used in constructing the feature."""
    
        def __new__(cls: Type["FeatureBaseCaller"], name: str, base: Any, attr: Any) -> "FeatureBaseCaller":
            """Initialize the metaclass.
    
            Args:
                name (str): The name of class.
                base (Any): Tuple of classes from which the current class derives.
                attr (Any): A dictionary that holds the namespaces for the class.
            Returns:
                FeatureBaseCaller : The metaclass.
            """
            if name != "FeatureBase":
                if "__call__" in attr:
                    wrapped = attr["__call__"]
                    attr["__call__"] = cls.call_decorate(wrapped)
                    functools.update_wrapper(attr["__call__"], wrapped, assigned=("__doc__", "__annotations__"), updated=())
            return super().__new__(cls, name, base, attr)
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseCaller.html#matlantis_features.features.base.FeatureBaseCaller.call_decorate)    @classmethod
        def call_decorate(cls, derived_call: Callable[..., T]) -> Callable[..., T]:
            """Decorate the __call__ function of the feature.
    
            Args:
                derived_call (... -> T): The original __call__ function of the feature.
            Returns:
                ... -> T : The decorated __call__ function.
            """
    
            def decorated(self: FeatureBase, *args: Any, **kwargs: Dict[str, Any]) -> T:
                args_list = [copy_to_child_conversion(arg) for arg in args]
                kwargs_new = {k: copy_to_child_conversion(v) for k, v in kwargs.items()}
                res = derived_call(self, *args_list, **kwargs_new)
                return res
    
            return decorated
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase)class FeatureBase(metaclass=FeatureBaseCaller):
        """The base of the matlantis features."""
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.__init__)    def __init__(self) -> None:
            """Initialize an instance."""
            super(FeatureBase, self).__init__()
            self._registered = False
            self._attached = False
            self._within_init_scope = False
            self._cost: Optional[FeatureCost] = None
            self._children: Set[str] = set()
            self._n_repeat = 1  # Currently this is always set to 1
            self._n_called = 0
            self._name = "root"
            self._total_name = self._name
    
            self._ctx: Optional[Context] = None
            self.run_id = 0
            self._tempdir = tempfile.TemporaryDirectory(prefix="matlantis_")
    
            try:
                s_id = os.getenv("GA4_SERVICE_ID")
                m_id = os.getenv("GA4_ID")
                c_id = os.getenv("MATLANTIS_USER_ID")
                send_data = os.getenv("MATLANTIS_COOKIE_CONSENT_PERFORMANCE") is not None
    
                if send_data and s_id is not None and m_id is not None and c_id is not None:
                    self.collector = GtagMP(s_id, m_id, c_id)
                    api_event = self.collector.create_new_event(name="matlantis_features_call")
                    api_event.set_event_param(name="feature_type", value=str(type(self).__name__))
                    self.collector.send(events=[api_event])
            except Exception:
                pass
    
    
    
        def __setattr__(self, name: str, value: Any) -> None:
            """Set attribution of the feature.
    
            Args:
                name (str): The name of attribution.
                value (Any): The value of the attribution.
            """
            if self.within_init_scope:
                if isinstance(value, FeatureBase):
                    if hasattr(self, name):
                        raise AttributeError(f"{name}: attribute exists")
                    value._name = name
                    value._total_name = self.total_name + "." + name
                    value.run_id = self.run_id
                    value._registered = True
                    self._children.add(name)
            else:
                if self.registered and not hasattr(self, name):
                    raise SetAttributeAfterRegisteredError
            super(FeatureBase, self).__setattr__(name, value)
    
        @property
        def name(self) -> str:
            """The name of the feature.
    
            Returns:
                str : The name of the feature.
            """
            return self._name
    
        @name.setter
        def name(self, name: str) -> None:
            """Set the name of the feature.
    
            Args:
                name (str): The name of the feature.
            """
            self._name = name
    
        @property
        def total_name(self) -> str:
            """The total name of the feature which is the combination of \
            own name and parent feature's name.
    
            Returns:
                str : The total name of the feature.
            """
            return self._total_name
    
        @property
        def n_repeat(self) -> int:
            """Get the maximum number of times that allowed to run the __call__ function.
    
            Returns:
                int : The maximum number of repeats.
            """
            return self._n_repeat
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.repeat)    def repeat(self: Self, n_repeat: int) -> Self:
            """Set the maximum number of times that allowed to run the __call__ function.
    
            Args:
                n_repeat (int): The maximum number of repeats.
            Returns:
                Self : The feature.
            """
            if self.registered:
                raise RepeatFeatureAfterRegsteredError
            assert self._n_repeat in (1, n_repeat)
            self._n_repeat = n_repeat
            return self
    
    
    
        @property
        def registered(self) -> bool:
            """Check whether the feature is registered or attached.
    
            Returns:
                bool : Whether the feature is registered or attached..
            """
            return bool(getattr(self, "_registered", False) or getattr(self, "_attached", False))
    
        @property
        def within_init_scope(self) -> bool:
            """Check whether or not in the init_scope.
    
            Returns:
                bool : Whether or not in the init_scope.
            """
            return bool(getattr(self, "_within_init_scope", False))
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.init_scope)    @contextlib.contextmanager
        def init_scope(self) -> Iterator[None]:
            """Context manager that enable to set attribution of the feature.
    
            Returns:
                Iterator[None] : Init_scope context manager.
            """
            old_flag = self.within_init_scope
            self._within_init_scope = True
            try:
                yield
            finally:
                self._within_init_scope = old_flag
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.attach_ctx)    def attach_ctx(self, ctx: Optional[Context] = None) -> None:
            """Attach the feature to matlantis_features.utils.Context.
    
            Args:
                ctx (Context or None, optional): The matlantis_features.utils.Context object. Defaults to None.
            """
            warnings.warn(
                "Context is not needed and is going to be deprecated in version 0.3.",
                FutureWarning,
            )
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.get_savedir_from_ctx)    def get_savedir_from_ctx(self) -> pathlib.Path:
            """Get the temporary save directory from the context.
    
            Returns:
                pathlib.Path : The temporary save directory .
            """
            return pathlib.Path(self._tempdir.name)
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.cost_estimate)    def cost_estimate(self, atoms: Optional[Union[ASEAtoms, MatlantisAtoms]] = None) -> FeatureCost:
            """Estimate the cost of the feature.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms or None, optional):
                  The input atoms. Defaults to None.
            Returns:
                FeatureCost : The cost of the feature.
            """
            if atoms is None:
                atoms = Atoms()
            total_cost: FeatureCost = self._internal_cost(len(atoms))
            for child_name in self._children:
                child = getattr(self, child_name)
                total_cost += child.cost_estimate(atoms) * child.n_repeat
            return total_cost
    
    
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            raise NotImplementedError
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the FeatureBase.
    
            Returns:
                dict[str, Any] : A dict containing a serialized FeatureBase.
            """
            raise NotImplementedError
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.check_estimator_fn)    def check_estimator_fn(self, estimator_fn: Optional[EstimatorFnType]) -> None:
            """
            Checks if the given estimator function is None and output a warning if so.
    
            Args:
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail. Defaults to None.
            """
            if estimator_fn is not None:
                return  # OK; estimator_fn is provided
    
            model_version_key = "MATLANTIS_PFP_MODEL_VERSION"
            calc_mode_key = "MATLANTIS_PFP_CALC_MODE"
            available_models = _default_calc_mode_table_by_model_version.keys()
    
            model_version = os.environ.get(model_version_key, Config.model_version)
            if model_version_key not in os.environ:
                warnings.warn(
                    f"Model version is not specified, default '{Config.model_version}' will be used."
                    f"Please specify it by providing an estimator_fn or setting the '{model_version_key}' environment variable.",
                    UserWarning,
                    stacklevel=2,
                )
            elif model_version not in available_models:
                version_info = get_version_info()
                pfp_api_client_version = version_info.get("pfp-api-client")
                raise ValueError(
                    f"Unexpected model version {model_version} was specified by '{model_version_key}' environment variable. "
                    f"Available model versions in pfp-api-client {pfp_api_client_version} are [{', '.join(available_models)}]."
                )
    
            calc_mode = os.environ.get(calc_mode_key)
            if calc_mode_key not in os.environ:
                calc_mode = _default_calc_mode_table_by_model_version[model_version]
                warnings.warn(
                    f"Calculation mode is not specified, default '{calc_mode}' will be used."
                    f"Please specify it by providing an estimator_fn or setting the '{calc_mode_key}' environment variable.",
                    UserWarning,
                    stacklevel=2,
                )
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBase.html#matlantis_features.features.base.FeatureBase.from_dict)    @classmethod
        def from_dict(cls, res: Dict[str, Any]) -> "FeatureBase":
            """Construct a FeatureBase object from serialized dict.
    
            Args:
                res (dict[str, Any]): A dict containing a serialized FeatureBase from to_dict().
            Returns:
                FeatureBase : The deserialized object from provided dict.
            """
            obj = from_dict(res)
            if not isinstance(obj, cls):
                raise TypeError(
                    f"Expected to build object of type {cls.__name__}, but input is a serialized object of type {obj.__class__.__name__}"
                )
            return obj
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput)@dataclass
    class FeatureBaseOutput:
        """The base of the matlantis features result output classes."""
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the FeatureBaseOutput.
    
            Returns:
                dict[str, Any] : A dict containing a serialized FeatureBaseOutput.
            """
            raw_dict = {f.name: getattr(self, f.name) for f in fields(self)}
            raw_dict["cls_path"] = f"{self.__class__.__module__}.{self.__class__.__name__}"
            return to_serialized_dict(raw_dict)
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseOutput.html#matlantis_features.features.base.FeatureBaseOutput.from_dict)    @classmethod
        def from_dict(cls, res: Dict[str, Any]) -> "FeatureBaseOutput":
            """Construct a FeatureBaseOutput object from serialized dict.
    
            Args:
                res (dict[str, Any]): A dict containing a serialized FeatureBaseOutput from to_dict().
            Returns:
                FeatureBaseOutput : The deserialized object from provided dict.
            """
            obj = from_dict(res)
            if not isinstance(obj, cls):
                raise TypeError(
                    f"Expected to build object of type {cls.__name__}, but input is a serialized object of type {obj.__class__.__name__}"
                )
            return obj
    
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult)@dataclass
    class FeatureBaseResult:
        """The base of the matlantis features result classes."""
    
        feature: FeatureBase
        call_params: Dict[str, Any]
        output: FeatureBaseOutput
        calculator_info: Optional[Dict[str, Any]]
        version_info: Dict[str, Any]
        user_info: Dict[str, Any]
    
        def __getattr__(self, name: str) -> Any:
            # Instead of calling result.output.force_constant,
            # simply call result.force_constant (to preserve backward compatibility).
            if name in [f.name for f in fields(self.output)]:
                obj = getattr(self.output, name)
                if isinstance(obj, dict):
                    return from_dict(self.output_state["name"])
                return obj
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult.to_dict)    def to_dict(self) -> Dict[str, Any]:
            """Dictionary representation of the FeatureBaseResult.
    
            Returns:
                dict[str, Any] : A dict containing a serialized FeatureBaseResult.
            """
            obj_dict = {
                "feature": self.feature,
                "call_params": self.call_params,
                "output": self.output,
                "calculator_info": self.calculator_info,
                "version_info": self.version_info,
                "user_info": self.user_info,
            }
            obj_dict["cls_path"] = f"{self.__class__.__module__}.{self.__class__.__name__}"
            return to_serialized_dict(obj_dict)
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult.from_dict)    @classmethod
        def from_dict(cls, res: Dict[str, Any]) -> "FeatureBaseResult":
            """Construct a FeatureBaseResult object from serialized dict.
    
            Args:
                res (dict[str, Any]): A dict containing a serialized FeatureBaseResult from to_dict().
            Returns:
                FeatureBaseResult : The deserialized object from provided dict.
            """
            obj = from_dict(res)
            if not isinstance(obj, cls):
                raise TypeError(
                    f"Expected to build object of type {cls.__name__}, but the input is a serialized object of type {obj.__class__.__name__}"
                )
            return obj
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult.save)    def save(self, filename: str) -> bool:
            """Construct a FeatureBaseResult object from serialized dict.
    
            Args:
                filename (str): Name of file to save the FeatureBaseResult to.
            Returns:
                bool : Whether saving the result was successful or not.
            """
            try:
                with open(filename, "w") as fd:
                    json.dump(self.to_dict(), fd, indent=4, cls=NpEncoder)
                    return True
            except Exception:
                return False
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.features.base.FeatureBaseResult.html#matlantis_features.features.base.FeatureBaseResult.load)    @classmethod
        def load(cls, filename: str) -> "FeatureBaseResult":
            """Construct a FeatureBaseResult object from serialized dict.
    
            Args:
                filename (str): Name of file to save the FeatureBaseResult to.
            Returns:
                FeatureBaseResult : The FeatureBaseResult object if loading was successful.
            """
            try:
                with open(filename, "r") as fd:
                    res = json.load(fd)
                if get_version_info() != res["version_info"]:
                    warnings.warn(
                        "Note: some package versions are different from those used at "
                        "the time when the save file was generated. This may cause "
                        "errors or other unexpected behavior when reloading.\n\n"
                        f"Current versions:\n{get_version_info()}\n\n"
                        f"Versions used at save file generation:\n{res['version_info']}"
                    )
                return cls.from_dict(res)
            except Exception as e:
                raise ValueError(f"Could not load result from {filename}: {e}")
    
    
    
    
    class RepeatFeatureAfterRegsteredError(Exception):
        """Change the setting of repeat of the feature after registration."""
    
        pass
    
    
    class RenameRegisteredFeatureError(Exception):
        """Change the name of the feature after registration."""
    
        pass
    
    
    class SetAttributeAfterRegisteredError(Exception):
        """Set the attribute of the feature after registration."""
    
        pass
    
    
    class NotAttachedContextError(Exception):
        """The feature is not attached to any matlantis_features.utils.Context object."""
    
        pass
    
    
    class CostEstimationWithoutContextError(Exception):
        """Estimate the cost of the feature when it is not attached to any \
        matlantis_features.utils.Context object."""
    
        pass
    
