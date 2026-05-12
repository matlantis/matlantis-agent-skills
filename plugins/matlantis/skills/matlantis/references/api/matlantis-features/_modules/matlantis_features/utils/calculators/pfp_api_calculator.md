# Source code for matlantis_features.utils.calculators.pfp_api_calculator
    
    
    import os
    from typing import Optional
    
    from ase.calculators.calculator import Calculator
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    from pfp_api_client.pfp.config import Config
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode, EstimatorMethodType
    
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    _calculator: Optional[ASECalculator] = None
    
    
    def _get_calc_mode(calc_mode: str) -> EstimatorCalcMode:
        calc_mode_obj = getattr(EstimatorCalcMode, calc_mode.upper(), None)
        if calc_mode_obj is None:
            raise ValueError(f"{calc_mode} is unknown calc mode.")
    
        return calc_mode_obj
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.utils.calculators.pfp_api_calculator.get_calculator.html#matlantis_features.utils.calculators.pfp_api_calculator.get_calculator)def get_calculator(
        estimator_fn: Optional[EstimatorFnType] = None,
    ) -> Calculator:
        """This function is not intended to be used by users."""
    
        if estimator_fn is not None:
            estimator = estimator_fn()
        else:
            model_version = os.environ.get("MATLANTIS_PFP_MODEL_VERSION", Config.model_version)
            calc_mode = os.environ.get("MATLANTIS_PFP_CALC_MODE", None)
            if calc_mode is not None:
                calc_mode = _get_calc_mode(calc_mode)
            estimator = Estimator(model_version=model_version, calc_mode=calc_mode)
    
        _calculator = ASECalculator(estimator)
        return _calculator
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.utils.calculators.pfp_api_calculator.pfp_estimator_fn.html#matlantis_features.utils.calculators.pfp_api_calculator.pfp_estimator_fn)def pfp_estimator_fn(
        *,
        max_retries: int = 10,
        model_version: str = Config.model_version,
        address: str = Config.address,
        priority: int = Config.priority,
        calc_mode: Optional[EstimatorCalcMode] = None,
        method_type: Optional[EstimatorMethodType] = None,
    ) -> EstimatorFnType:
        """An estimator function for a PFP estimator.
    
        For more information about the estimator function,
        please read our `guidebook(en) </api/resource/documents/matlantis-guidebook/en/
        about-matlantis-features.html#customizing-estimator-used-in-matlantis-features>`_
        or `guidebook(ja) </api/resource/documents/matlantis-guidebook/ja/
        about-matlantis-features.html#matlantis-featuresestimator>`_.
    
        Args:
            max_retries (int, optional):
                Maximum number of retries when pfp-api-client fails to
                create an estimator. Defaults to 10.
            model_version (str, optional):
                The model version to be set.
            address (str, optional):
                The address of PFP API server.
            priority (int, optional):
                An integer (1 to 100) that determines how to prioritize requests
                when there is a lack of Token Availability. Defaults to 100.
            secure (bool, optional):
                The flag to control if it uses SSL.
            ca_cert_path (str or None, optional):
                The path to certificate file
            calc_mode (EstimatorCalcMode or None, optional):
                The calculation mode.
            method_type (EstimatorMethodType or None, optional):
                The inference method type to be set.
        Returns:
            EstimatorFnType : callable object that returns an estimator.
        """
    
        def estimator_fn() -> Estimator:
            return Estimator(
                max_retries=max_retries,
                model_version=model_version,
                address=address,
                priority=priority,
                calc_mode=calc_mode,
                method_type=method_type,
            )
    
        return estimator_fn
    
    
    
