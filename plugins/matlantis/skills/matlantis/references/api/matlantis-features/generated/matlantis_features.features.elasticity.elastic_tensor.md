# matlantis_features.features.elasticity.elastic_tensor#

Features

[`ElasticTensorFeature`](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature "matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature")(diagonal_strains, ...) | The matlantis-feature for calculating the elastic tensor for the specified atoms, by applying the specified strain.  
---|---  
  
Classes

[`ElasticTensor`](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")(stress_data, strain_data) | Elastic tensor object calculated from the stress-strain relation.  
---|---  
  
Functions

[`make_deformed_atoms`](matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms.html#matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms "matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms")(atoms, dmat) | Create deformed Atoms object by deformation tensor.  
---|---  
[`stress_from_deform`](matlantis_features.features.elasticity.elastic_tensor.stress_from_deform.html#matlantis_features.features.elasticity.elastic_tensor.stress_from_deform "matlantis_features.features.elasticity.elastic_tensor.stress_from_deform")(dmat, atoms, opt_fn, ...) | Calculate the stress tensor from the target atoms and deformation tensor, using specified optimization method.  
  
[ __ previous matlantis_features.features.elasticity.elastic_properties.expand_units ](matlantis_features.features.elasticity.elastic_properties.expand_units.html "previous page") [ next matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature __](matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html "next page")
