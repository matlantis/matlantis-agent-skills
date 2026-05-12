# matlantis_features.features.common.opt#

Features

[`ASEOptFeature`](matlantis_features.features.common.opt.ASEOptFeature.html#matlantis_features.features.common.opt.ASEOptFeature "matlantis_features.features.common.opt.ASEOptFeature")(get_optimizer[, n_run, fmax, ...]) | The matlantis-feature that wraps around an ASE Optimizer object.  
---|---  
[`BFGSASEOptFeature`](matlantis_features.features.common.opt.BFGSASEOptFeature.html#matlantis_features.features.common.opt.BFGSASEOptFeature "matlantis_features.features.common.opt.BFGSASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the ASE implementation of the BFGS algorithm.  
[`BFGSLineSearchASEOptFeature`](matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature "matlantis_features.features.common.opt.BFGSLineSearchASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the ASE implementation of the BFGS line search algorithm.  
[`FireASEOptFeature`](matlantis_features.features.common.opt.FireASEOptFeature.html#matlantis_features.features.common.opt.FireASEOptFeature "matlantis_features.features.common.opt.FireASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the Fast Inertial Relaxation Engine (FIRE), as implemented within ASE.  
[`FireLBFGSASEOptFeature`](matlantis_features.features.common.opt.FireLBFGSASEOptFeature.html#matlantis_features.features.common.opt.FireLBFGSASEOptFeature "matlantis_features.features.common.opt.FireLBFGSASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that combines FIRE and LBFGS optimization algorithms.  
[`LBFGSASEOptFeature`](matlantis_features.features.common.opt.LBFGSASEOptFeature.html#matlantis_features.features.common.opt.LBFGSASEOptFeature "matlantis_features.features.common.opt.LBFGSASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the ASE implementation of the LBFGS algorithm.  
[`LBFGSLineSearchASEOptFeature`](matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature.html#matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature "matlantis_features.features.common.opt.LBFGSLineSearchASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the ASE implementation of the BFGS line search algorithm.  
[`MDMinASEOptFeature`](matlantis_features.features.common.opt.MDMinASEOptFeature.html#matlantis_features.features.common.opt.MDMinASEOptFeature "matlantis_features.features.common.opt.MDMinASEOptFeature")([n_run, fmax, ...]) | The matlantis-feature that uses the ASE implementation of the MDMin algorithm.  
  
FeatureResults

[`OptFeatureResult`](matlantis_features.features.common.opt.OptFeatureResult.html#matlantis_features.features.common.opt.OptFeatureResult "matlantis_features.features.common.opt.OptFeatureResult")(feature, call_params, ...) | Class with the results for various optimizer features.  
---|---  
  
Classes

[`OptFeatureBase`](matlantis_features.features.common.opt.OptFeatureBase.html#matlantis_features.features.common.opt.OptFeatureBase "matlantis_features.features.common.opt.OptFeatureBase")() | Template base class for optimizer features.  
---|---  
[`OptFeatureOutput`](matlantis_features.features.common.opt.OptFeatureOutput.html#matlantis_features.features.common.opt.OptFeatureOutput "matlantis_features.features.common.opt.OptFeatureOutput")(atoms, converged, steps, ...) | Class containing the outputs for results of various optimizer features.  
  
[ __ previous OPT Features ](../matlantis_features.features.common.opt.html "previous page") [ next matlantis_features.features.common.opt.ASEOptFeature __](matlantis_features.features.common.opt.ASEOptFeature.html "next page")
