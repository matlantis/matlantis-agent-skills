# PFP#

Module containing main PFP functionalities.

## Estimator#

This module contains the core estimator class for inference using the PFP API.

[`Estimator`](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator "pfp_api_client.pfp.estimator.Estimator")([max_retries, model_version, ...]) | Core class of PFP.  
---|---  
[`EstimatorCalcMode`](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorCalcMode "pfp_api_client.pfp.estimator.EstimatorCalcMode")(value[, names, module, ...]) | Enum class which is used to determine calc_mode in Estimator.  
[`EstimatorMethodType`](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorMethodType "pfp_api_client.pfp.estimator.EstimatorMethodType")(value[, names, module, ...]) | Enum class which is used to determine method_type in Estimator.  
[`EstimatorSystem`](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.EstimatorSystem "pfp_api_client.pfp.estimator.EstimatorSystem")(properties, coordinates, cell) | Input argument of Estimator.estimate function.  
  
## Calculator#

This module allows the PFP estimator to be used as a calculator within ASE.

[`ASECalculator`](pfp_api_client.pfp.calculators.html#pfp_api_client.pfp.calculators.ase_calculator.ASECalculator "pfp_api_client.pfp.calculators.ase_calculator.ASECalculator")(estimator) | ASECalculator is a derivation of ase.calculators.calculator.  
---|---  
  
[ __ previous PFP API Client Reference ](pfp_api_client.html "previous page") [ next Utilities __](pfp_api_client.pfp.utils.html "next page")

__On this page

  * Estimator 
  * Calculator 


