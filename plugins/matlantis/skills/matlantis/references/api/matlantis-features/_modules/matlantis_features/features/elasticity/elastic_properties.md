# Source code for matlantis_features.features.elasticity.elastic_properties
    
    
    import itertools
    from typing import Dict, List, Optional, Union
    
    import numpy as np
    from ase import units
    
    from ...utils import FeatureCost
    from ..base import FeatureBase
    from .elastic_tensor import ElasticTensor
    from .utils import ij_to_voigt, to_unit_vec
    
    
    def _s_voigt_coeff(p: int, q: int) -> float:
        return 1.0 / ((1 + p // 3) * (1 + q // 3))
    
    
    def _calc_smat(compl_tensor: np.ndarray) -> np.ndarray:
        # Calculate "Smat" from compliance tensor. See:
        #   https://github.com/coudertlab/elate
        #   https://github.com/coudertlab/elate/blob/master/elastic.py#L624
        smat = np.zeros((3, 3, 3, 3))
        for i, j, k, l in itertools.product(range(3), range(3), range(3), range(3)):
            smat[i, j, k, l] = _s_voigt_coeff(ij_to_voigt(i, j), ij_to_voigt(k, l)) * compl_tensor[ij_to_voigt(i, j), ij_to_voigt(k, l)]
        return smat  # type: ignore
    
    
    PressureResult = Dict[str, float]
    PostElasticFeatureResult = Dict[str, Union[float, PressureResult, Dict[str, List[float]]]]
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.expand_units.html#matlantis_features.features.elasticity.elastic_properties.expand_units)def expand_units(atom_unit_pressure: float) -> PressureResult:
        """Expand pressure units.
    
        Args:
            atom_unit_pressure (float): Input pressure results in GPa unit.
        Returns:
            PressureResult : Expanded result.
        """
        # Default unit for pressure is GPa ?
        return {
            "eV/A^3": atom_unit_pressure * units.GPa,
            "GPa": atom_unit_pressure,
            "bar": atom_unit_pressure / units.GPa * units.bar,
        }
    
    
    
    
    def _dir_to_str(dir: List[float]) -> str:
        return f"[{dir[0]} {dir[1]} {dir[2]}]"
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties)class ElasticProperties(object):
        """Elastic properties."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.__init__)    def __init__(self, data: np.ndarray):
            """Initialize an instance.
    
            Args:
                data (np.ndarray): The elastic tensor.
            """
            if data.shape != (6, 6):
                raise ValueError(f"invalid elastic_tensor shape: {data.shape}")
            self._elastic_tensor = data
            self._compliance_tensor = np.linalg.inv(self.elastic_tensor)
    
            self._Smat = _calc_smat(self._compliance_tensor)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.from_elastic_tensor)    @classmethod
        def from_elastic_tensor(cls, elastic_tensor: ElasticTensor) -> "ElasticProperties":
            """Create elastic properties object from elastic tensor.
    
            Args:
                elastic_tensor (ElasticTensor): Elastic tensor.
            Returns:
                ElasticProperties : Resulting object.
            """
            return cls(elastic_tensor.elastic_tensor)
    
    
    
        @property
        def elastic_tensor(self) -> np.ndarray:
            """Elastic tensor.
    
            Returns:
                np.ndarray : Numpy ndarray representing elastic tensor.
            """
            return self._elastic_tensor
    
        @property
        def Smat(self) -> np.ndarray:
            """Smat.
    
            Returns:
                np.ndarray : Numpy ndarray representing Smat.
            """
            return self._Smat
    
        @property
        def compliance_tensor(self) -> np.ndarray:
            """Calculate compliance tensor, s.
    
            Returns:
                np.ndarray : Numpy ndarray representing compliance tensor.
            """
            return self._compliance_tensor  # type: ignore
    
        @property
        def bulk_modulus_voigt(self) -> float:
            """Calculate Voigt bulk modulus, K_v.
    
            9 * K_v = (C_11 + C_22 + C_33) + 2 * (C_12 + C_23 + C_13)
    
            Returns:
                float : bulk modulus value.
            """
            return float(np.mean(self.elastic_tensor[:3, :3]))
    
        @property
        def bulk_modulus_reuss(self) -> float:
            """Calculate Reuss bulk modulus, K_r.
    
            1/K_r = (s_11 + s_22 + s_33) + 2 * (s_12 + s_23 + s_13)
    
            Returns:
                float : bulk modulus value.
            """
            return 1.0 / float(np.sum(self.compliance_tensor[:3, :3]))
    
        @property
        def bulk_modulus(self) -> float:
            """Calculate Voigt-Reuss-Hill bulk modulus, K_vrh.
    
            k_vrh = 0.5 * (K_v + K_r)
    
            Returns:
                float : bulk modulus value.
            """
            return 0.5 * (self.bulk_modulus_voigt + self.bulk_modulus_reuss)
    
        @property
        def shear_modulus_voigt(self) -> float:
            """Calculate Voigt shear modulus, G_v.
    
            15 * G_v = (C_11 + C_22 + C_33) - (C_12 + C_23 + C_13) + 3 * (C_44 + C_55 + C_66)
    
            Returns:
                float : shear modulus value.
            """
            return (
                2.0 * float(self.elastic_tensor[:3, :3].trace())
                - float(np.triu(self.elastic_tensor[:3, :3]).sum())
                + 3 * float(self.elastic_tensor[3:, 3:].trace())
            ) / 15.0
    
        @property
        def shear_modulus_reuss(self) -> float:
            """Calculate Reuss shear modulus, G_r.
    
            15 / G_v = 4 * (s_11 + s_22 + s_33) - 4 * (s_12 + s_23 + s_13) + 3 * (s_44 + s_55 + s_66)
    
            Returns:
                float : shear modulus value.
            """
            return 15.0 / (
                8.0 * float(self.compliance_tensor[:3, :3].trace())
                - 4.0 * float(np.triu(self.compliance_tensor[:3, :3]).sum())
                + 3.0 * float(self.compliance_tensor[3:, 3:].trace())
            )
    
        @property
        def shear_modulus(self) -> float:
            """Calculate Voigt-Reuss-Hill shear modulus, G_vrh.
    
            G_vrh = 0.5 * (G_v + G_r)
    
            Returns:
                float : shear modulus value.
            """
            return 0.5 * (self.shear_modulus_voigt + self.shear_modulus_reuss)
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.get_directional_shear_modulus)    def get_directional_shear_modulus(self, vec1: np.ndarray, vec2: np.ndarray) -> float:
            """Calculate directional shear modulus.
    
            Adopted the method from elate https://github.com/coudertlab/elate
            shear modulus along specific directions.
    
            Args:
                vec1 (np.ndarray): The first 3D vector specifying the direction.
                vec2 (np.ndarray): The second 3D vector specifying the direction.
                                   vec2 should be perpendicular to vec1.
            Returns:
                float : Resulting directional shear modulus.
            """
            vec1 = to_unit_vec(vec1)
            vec2 = to_unit_vec(vec2)
            if not np.abs(np.dot(vec1, vec2)) < 1e-3:
                raise ValueError("n and m must be orthogonal")
    
            r = float(np.einsum("i,j,k,l,ijkl", vec1, vec2, vec1, vec2, self.Smat))
            return 1 / (4 * r)
    
    
    
        @property
        def poisson_ratio(self) -> float:
            """Calculate Poisson's ratio, mu.
    
            (3 * K_vrh + 2 * G_vrh) / (6 * K_vrh + 2 * G_vrh)
    
            Returns:
                float : Poisson's ratio value.
            """
            return (3 * self.bulk_modulus - 2 * self.shear_modulus) / (6 * self.bulk_modulus + 2 * self.shear_modulus)
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.get_directional_poisson_ratio)    def get_directional_poisson_ratio(self, vec1: np.ndarray, vec2: np.ndarray) -> float:
            """Calculate directional Poisson's ratio.
    
            Adopted the method from elate https://github.com/coudertlab/elate
            poisson ratio along specific directions.
    
            Args:
                vec1 (np.ndarray): First 3D vector specifying the direction.
                vec2 (np.ndarray): Second 3D vector specifying the direction.
                                   vec2 should be perpendicular to vec1.
            Returns:
                float : Resulting directional Poisson's ratio.
            """
            vec1 = to_unit_vec(vec1)
            vec2 = to_unit_vec(vec2)
            if not np.abs(np.dot(vec1, vec2)) < 1e-3:
                raise ValueError("n and m must be orthogonal")
    
            r1 = float(np.einsum("i,j,k,l,ijkl", vec1, vec1, vec2, vec2, self.Smat))
            r2 = float(np.einsum("i,j,k,l,ijkl", vec1, vec1, vec1, vec1, self.Smat))
            return -1.0 * r1 / r2
    
    
    
        @property
        def anisotropy(self) -> float:
            """Calculate universal anisotropy.
    
            5 * (G_v/G_r) + (K_v/K_r) - 6
            Proposed by S. I. Ranganathan and M. Ostoja-Starzewski
            (Phys. Rev. Lett. 101, 3, 2008).
    
            Returns:
                float : Resulting universal anisotropy value.
            """
            return 5.0 * self.shear_modulus_voigt / self.shear_modulus_reuss + self.bulk_modulus_voigt / self.bulk_modulus_reuss - 6.0
    
        @property
        def young_modulus(self) -> float:
            """Calculate Young's modulus.
    
            9 * K_vrh * G_vrh / (3 * K_vrh + G_vrh)
    
            Returns:
                float : Resulting Young's modulus value.
            """
            return 9 * self.bulk_modulus * self.shear_modulus / (3 * self.bulk_modulus + self.shear_modulus)
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.get_directional_young_modulus)    def get_directional_young_modulus(self, vec: np.ndarray) -> float:
            """Calculate directional Young's modulus. \
            Adopted the method from elate https://github.com/coudertlab/elate \
            young's modulus along specific direction.
    
            Args:
                vec (np.ndarray): 3D vector specifying the direction.
    
            Returns:
                float : Resulting directional Young's modulus value.
            """
            vec = to_unit_vec(vec)
    
            r = float(np.einsum("i,j,k,l,ijkl", vec, vec, vec, vec, self.Smat))
            return 1 / r
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.get_directional_linear_compressibility)    def get_directional_linear_compressibility(self, vec: np.ndarray) -> float:
            """Calculate directional linear compressibility. \
            Adopted the method from elate https://github.com/coudertlab/elate \
            linear compressibility along specific direction.
    
            Args:
                vec (np.ndarray): 3D vector specifying the direction.
    
            Returns:
                float : Resulting directional linear compressibility value.
            """
            vec = to_unit_vec(vec)
    
            r = float(np.einsum("i,j,ijkk", vec, vec, self.Smat))
            return 1000 * r
    
    
    
        @property
        def vicker_hardness(self) -> float:
            """Calculate vicker hardness H_v.
    
            H_v = 0.92 * (G/K)^1.137 * G^0.708
            (Liu et al. Phys. Rev. B 90, 134102, 2014)
    
            Returns:
                float : Resulting vicker hardness value.
            """
            return (  # type: ignore[no-any-return]
                0.92 * (self.shear_modulus / self.bulk_modulus) ** 1.137 * self.shear_modulus**0.708
            )
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.ElasticProperties.html#matlantis_features.features.elasticity.elastic_properties.ElasticProperties.get_prop_results)    def get_prop_results(self, dir1: Optional[List[float]] = None, dir2: Optional[List[float]] = None) -> PostElasticFeatureResult:
            """Calculate elastic properties.
    
            Args:
                dir1 (list[float] or None, optional): First 3D vector specifying the direction.
                    If None, directional properties will not be calculated.
                    Defaults to None.
                dir2 (list[float] or None, optional): Second 3D vector specifying the direction.
                    vec2 should be perpendicular to vec1.
                    Defaults to None.
            Returns:
                PostElasticFeatureResult : Result object containing the calculated properties.
            """
            results: PostElasticFeatureResult = {
                "poisson_ratio": self.poisson_ratio,
                "anisotropy": self.anisotropy,
                "bulk_modulus": expand_units(self.bulk_modulus),
                "young's_modulus": expand_units(self.young_modulus),
                "shear_modulus": expand_units(self.shear_modulus),
            }
            return_elems = [
                [0, 0],
                [1, 1],
                [2, 2],
                [0, 1],
                [0, 2],
                [1, 2],
                [3, 3],
                [4, 4],
                [5, 5],
            ]
            for i, j in return_elems:
                results[f"C{i + 1}{j + 1}"] = expand_units(float(self.elastic_tensor[i, j]))
    
            if dir1 is not None and dir2 is not None:
                dir1_str = _dir_to_str(dir1)
                dir2_str = _dir_to_str(dir2)
                results["config"] = {
                    "direction": dir1,
                    "direction2": dir2,
                }
    
                # TODO(teaism_compatibility) reasonable keyname
                results[f"young's modulus along {dir1_str}"] = expand_units(self.get_directional_young_modulus(np.array(dir1)))
    
                results[f"linear compressibility along {dir1_str}"] = self.get_directional_linear_compressibility(np.array(dir1))
    
                results[f"shear modulus along {dir1_str} {dir2_str}"] = expand_units(
                    self.get_directional_shear_modulus(np.array(dir1), np.array(dir2))
                )
                results[f"poisson ratio along {dir1_str} {dir2_str}"] = self.get_directional_poisson_ratio(np.array(dir1), np.array(dir2))
    
            return results
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature.html#matlantis_features.features.elasticity.elastic_properties.PostElasticPropertiesFeature)class PostElasticPropertiesFeature(FeatureBase):
        """The matlantis-feature for calculating elastic properties.
    
        Including bulk modulus, \
        shear modulus, Young's modulus, poisson ratio etc.
        """
    
        def __call__(self, elastic_tensor: ElasticTensor) -> ElasticProperties:
            """Call function.
    
            Args:
                elastic_tensor (ElasticTensor): Elastic tensor object calculated by ElasticTensorFeature.
            Returns:
                ElasticProperties : Resulting elastic properties object.
            """
            props = ElasticProperties.from_elastic_tensor(elastic_tensor)
            return props
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=10.0)
    
    
    
