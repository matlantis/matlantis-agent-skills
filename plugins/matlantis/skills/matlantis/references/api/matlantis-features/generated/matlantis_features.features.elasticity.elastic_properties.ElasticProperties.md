# matlantis_features.features.elasticity.elastic_properties.ElasticProperties#

_class _matlantis_features.features.elasticity.elastic_properties.ElasticProperties(_data : ndarray_)[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties)#
    

Bases: `object`

Elastic properties.

Methods

`__init__`(data) | Initialize an instance.  
---|---  
`from_elastic_tensor`(elastic_tensor) | Create elastic properties object from elastic tensor.  
`get_directional_linear_compressibility`(vec) | Calculate directional linear compressibility.  
`get_directional_poisson_ratio`(vec1, vec2) | Calculate directional Poisson's ratio.  
`get_directional_shear_modulus`(vec1, vec2) | Calculate directional shear modulus.  
`get_directional_young_modulus`(vec) | Calculate directional Young's modulus.  
`get_prop_results`([dir1, dir2]) | Calculate elastic properties.  
  
__init__(_data : ndarray_)[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.__init__)#
    

Initialize an instance.

Parameters
    

**data** (_np.ndarray_) – The elastic tensor.

_classmethod _from_elastic_tensor(_elastic_tensor : [ElasticTensor](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")_) → ElasticProperties[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.from_elastic_tensor)#
    

Create elastic properties object from elastic tensor.

Parameters
    

**elastic_tensor** ([_ElasticTensor_](matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor "matlantis_features.features.elasticity.elastic_tensor.ElasticTensor")) – Elastic tensor.

Returns
    

Resulting object.

Return type
    

ElasticProperties

get_directional_linear_compressibility(_vec : ndarray_) → float[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.get_directional_linear_compressibility)#
    

Calculate directional linear compressibility. Adopted the method from elate <https://github.com/coudertlab/elate> linear compressibility along specific direction.

Parameters
    

**vec** (_np.ndarray_) – 3D vector specifying the direction.

Returns
    

Resulting directional linear compressibility value.

Return type
    

float

get_directional_poisson_ratio(_vec1 : ndarray_, _vec2 : ndarray_) → float[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.get_directional_poisson_ratio)#
    

Calculate directional Poisson’s ratio.

Adopted the method from elate <https://github.com/coudertlab/elate> poisson ratio along specific directions.

Parameters
    

  * **vec1** (_np.ndarray_) – First 3D vector specifying the direction.

  * **vec2** (_np.ndarray_) – Second 3D vector specifying the direction. vec2 should be perpendicular to vec1.



Returns
    

Resulting directional Poisson’s ratio.

Return type
    

float

get_directional_shear_modulus(_vec1 : ndarray_, _vec2 : ndarray_) → float[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.get_directional_shear_modulus)#
    

Calculate directional shear modulus.

Adopted the method from elate <https://github.com/coudertlab/elate> shear modulus along specific directions.

Parameters
    

  * **vec1** (_np.ndarray_) – The first 3D vector specifying the direction.

  * **vec2** (_np.ndarray_) – The second 3D vector specifying the direction. vec2 should be perpendicular to vec1.



Returns
    

Resulting directional shear modulus.

Return type
    

float

get_directional_young_modulus(_vec : ndarray_) → float[[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.get_directional_young_modulus)#
    

Calculate directional Young’s modulus. Adopted the method from elate <https://github.com/coudertlab/elate> young’s modulus along specific direction.

Parameters
    

**vec** (_np.ndarray_) – 3D vector specifying the direction.

Returns
    

Resulting directional Young’s modulus value.

Return type
    

float

get_prop_results(_dir1 : Optional[List[float]] = None_, _dir2 : Optional[List[float]] = None_) → Dict[str, Union[float, Dict[str, float], Dict[str, List[float]]]][[source]](../_modules/matlantis_features/features/elasticity/elastic_properties.html#ElasticProperties.get_prop_results)#
    

Calculate elastic properties.

Parameters
    

  * **dir1** (_list_ _[__float_ _] or_ _None_ _,__optional_) – First 3D vector specifying the direction. If None, directional properties will not be calculated. Defaults to None.

  * **dir2** (_list_ _[__float_ _] or_ _None_ _,__optional_) – Second 3D vector specifying the direction. vec2 should be perpendicular to vec1. Defaults to None.



Returns
    

Result object containing the calculated properties.

Return type
    

PostElasticFeatureResult

[ __ previous matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature ](matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature.html "previous page") [ next matlantis_features.features.elasticity.elastic_properties.expand_units __](matlantis_features.features.elasticity.elastic_properties.expand_units.html "next page")
