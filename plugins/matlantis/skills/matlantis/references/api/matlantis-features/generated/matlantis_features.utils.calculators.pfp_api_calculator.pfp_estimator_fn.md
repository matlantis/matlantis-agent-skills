# matlantis_features.utils.calculators.pfp_api_calculator.pfp_estimator_fn#

matlantis_features.utils.calculators.pfp_api_calculator.pfp_estimator_fn(_*_ , _max_retries : int = 10_, _model_version : str = 'v8.0.0'_, _address : str = 'localhost:5000'_, _priority : int = 100_, _calc_mode : Optional[EstimatorCalcMode] = None_, _method_type : Optional[EstimatorMethodType] = None_) → Callable[[], Estimator][[source]](../_modules/matlantis_features/utils/calculators/pfp_api_calculator.html#pfp_estimator_fn)#
    

An estimator function for a PFP estimator.

For more information about the estimator function, please read our [guidebook(en)](/api/resource/documents/matlantis-guidebook/en/about-matlantis-features.html#customizing-estimator-used-in-matlantis-features) or [guidebook(ja)](/api/resource/documents/matlantis-guidebook/ja/about-matlantis-features.html#matlantis-featuresestimator).

Parameters
    

  * **max_retries** (_int_ _,__optional_) – Maximum number of retries when pfp-api-client fails to create an estimator. Defaults to 10.

  * **model_version** (_str_ _,__optional_) – The model version to be set.

  * **address** (_str_ _,__optional_) – The address of PFP API server.

  * **priority** (_int_ _,__optional_) – An integer (1 to 100) that determines how to prioritize requests when there is a lack of Token Availability. Defaults to 100.

  * **secure** (_bool_ _,__optional_) – The flag to control if it uses SSL.

  * **ca_cert_path** (_str_ _or_ _None_ _,__optional_) – The path to certificate file

  * **calc_mode** (_EstimatorCalcMode_ _or_ _None_ _,__optional_) – The calculation mode.

  * **method_type** (_EstimatorMethodType_ _or_ _None_ _,__optional_) – The inference method type to be set.



Returns
    

callable object that returns an estimator.

Return type
    

EstimatorFnType

[ __ previous matlantis_features.utils.calculators.pfp_api_calculator.get_calculator ](matlantis_features.utils.calculators.pfp_api_calculator.get_calculator.html "previous page") [ next matlantis_features.functions __](../matlantis_features.functions.html "next page")
