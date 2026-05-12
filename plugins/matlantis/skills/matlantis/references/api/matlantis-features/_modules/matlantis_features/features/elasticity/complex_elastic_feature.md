# Source code for matlantis_features.features.elasticity.complex_elastic_feature
    
    
    import logging
    from dataclasses import dataclass
    from typing import Union
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.common.opt import ASEOptFeature, OptFeatureResult
    from matlantis_features.features.elasticity import PostElasticPropertiesFeature
    from matlantis_features.features.elasticity.elastic_properties import ElasticProperties
    from matlantis_features.features.elasticity.elastic_tensor import (
        ElasticTensor,
        ElasticTensorFeature,
    )
    from matlantis_features.utils import FeatureCost
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult.html#matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeatureResult)@dataclass
    class ComplexElasticFeatureResult:
        """A dataclass for result of ComplexElasticFeature."""
    
        run_opt_feature_result: OptFeatureResult
        elastic_tensor: ElasticTensor
        elastic_properties: ElasticProperties
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature.html#matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature)class ComplexElasticFeature(FeatureBase):
        """The matlantis-feature for calculating elasticity related properties.
    
        This feature calls both ElasticTensorFeature and PostElasticPropertiesFeature.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature.html#matlantis_features.features.elasticity.complex_elastic_feature.ComplexElasticFeature.__init__)    def __init__(
            self,
            run_ase_opt_feature: ASEOptFeature,
            run_elastic_tensor_feature: ElasticTensorFeature,
            post_elastic_feature: PostElasticPropertiesFeature,
        ) -> None:
            """Initialize an instance.
    
            Args:
                run_ase_opt_feature (ASEOptFeature):
                  Optimizer feature object used for elastic tensor calculation.
                run_elastic_tensor_feature (ElasticTensorFeature):
                  Feature object for calculating the elastic tensor.
                post_elastic_feature (PostElasticPropertiesFeature):
                  Feature object for calculating the elastic properties.
            """
            super(ComplexElasticFeature, self).__init__()
            with self.init_scope():
                self.run_ase_opt_feature = run_ase_opt_feature
                self.run_elastic_tensor_feature = run_elastic_tensor_feature
                self.post_elastic_feature = post_elastic_feature
    
                if self.run_ase_opt_feature.filter is None:
                    logging.warning(
                        "Cell optimization is usually needed to get correct elasticity property. "
                        "The cell optimization can be enabled by set the 'filter=True' in the "
                        "ASEOPTFeature that passed to the 'run_ase_opt_feature'."
                    )
    
    
    
        def __call__(self, atoms: Union[ASEAtoms, MatlantisAtoms]) -> ComplexElasticFeatureResult:
            """Perform the feature's calculation.
    
            Args:
                atoms (ASEAtoms or MatlantisAtoms):
                  Atoms object containing the target system.
            Returns:
                ComplexElasticFeatureResult:
                  Dataclass object containing the calculation result.
            """
            run_opt_feature_result = self.run_ase_opt_feature(atoms)
            elastic_tensor = self.run_elastic_tensor_feature(run_opt_feature_result.atoms)
            elastic_properties = self.post_elastic_feature(elastic_tensor)
    
            return ComplexElasticFeatureResult(
                run_opt_feature_result=run_opt_feature_result,
                elastic_tensor=elastic_tensor,
                elastic_properties=elastic_properties,
            )
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
