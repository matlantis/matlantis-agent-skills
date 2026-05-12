# Source code for matlantis_features.features.elasticity.tensors
    
    
    from typing import List
    
    import numpy as np
    import scipy
    
    from .utils import ij_to_voigt, is_voigt_diag, map2d, select, swap_ij, voigt_to_ij
    
    
    class _TensorVoigt:
        def __init__(self, value: np.ndarray) -> None:
            """Tensor in Voigt's notation.
    
            Args:
                value (np.ndarray):
                  numpy array of (6,) shape containing the tensor
                  elements in the Voigt's notation.
            """
            self.voigt = True
            if value.shape != (6,):
                raise ValueError(f"invalid value shape {value.shape}")
            self._data_voigt = value
    
        def set_by_voigt(self, vix: int, value: float) -> None:
            """Set tensor element by voigt index.
    
            Args:
                vix (int): Voigt index (0-5).
                value (float): Value to be stored in the tensor.
            """
            self._data_voigt[vix] = value
    
        def get_by_voigt(self, vix: int) -> float:
            """Get tensor element by voigt index.
    
            Args:
                vix (int): Voigt index (0-5).
            Returns:
                float : Value of the tensor of the specified index.
            """
            result: float = self._data_voigt[vix]
            return result
    
        def to_3x3_ndarray(self) -> np.ndarray:
            """Get 3x3 matrix representation of this tensor.
    
            Returns:
                np.ndarray : numpy array of (3,3) shape.
            """
            result = np.empty(shape=(3, 3))
            for i in range(3):
                for j in range(3):
                    result[i, j] = self._data_voigt[ij_to_voigt(i, j)]
            return result  # type: ignore
    
        def to_voigt_ndarray(self) -> np.ndarray:
            """Get 6x1 Voigt representation of this tensor.
    
            Returns:
                np.ndarray : numpy array of (6,1) shape.
            """
            return self._data_voigt
    
        def __repr__(self) -> str:
            """Get string representation.
    
            Returns:
                str : string representation of this tensor.
            """
            return f"_TensorVoigt(\n{self._data_voigt}\n)"
    
    
    class _Tensor3x3:
        def __init__(self, value: np.ndarray) -> None:
            """Tensor in 3x3 symmetric matrix representation.
    
            Args:
                value (np.ndarray):
                  numpy array of (3,3) shape containing the tensor
                  elements in the symmetrix matrix representation.
            """
            self.voigt = False
            if value.shape != (3, 3):
                raise ValueError(f"invalid value shape {value.shape}")
            self._data_3x3 = value
    
        def set_by_voigt(self, vix: int, value: float) -> None:
            """Set tensor element by voigt index.
    
            Args:
                vix (int): Voigt index (0-5).
                value (float): Value to be stored in the tensor.
            """
            i_j = voigt_to_ij(vix)
            self._data_3x3[i_j] = value
            if not is_voigt_diag(vix):
                self._data_3x3[swap_ij(i_j)] = value
    
        def get_by_voigt(self, vix: int) -> float:
            """Get tensor element by voigt index.
    
            Args:
                vix (int): Voigt index (0-5).
            Returns:
                float : Value of the tensor of the specified index.
            """
            result: float = self._data_3x3[voigt_to_ij(vix)]
            return result
    
        def to_3x3_ndarray(self) -> np.ndarray:
            """Get 3x3 matrix representation of this tensor.
    
            Returns:
                np.ndarray : numpy array of (3,3) shape.
            """
            return self._data_3x3
    
        def to_voigt_ndarray(self) -> np.ndarray:
            """Get 6x1 Voigt representation of this tensor.
    
            Returns:
                np.ndarray : numpy array of (6,1) shape.
            """
            result = np.empty(shape=(6,))
            for vix in range(6):
                result[vix] = self._data_3x3[voigt_to_ij(vix)]
            return result  # type: ignore
    
        def __repr__(self) -> str:
            """Get string representation.
    
            Returns:
                str : string representation of this tensor.
            """
            return f"_Tensor3x3(\n{self._data_3x3}\n)"
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor)class StrainTensor(_Tensor3x3):
        """Strain tensor object."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensor.html#matlantis_features.features.elasticity.tensors.StrainTensor.__init__)    def __init__(self, vix: int, strain: float) -> None:
            """Initialize an instance.
    
            Args:
                vix (int): Voigt index.
                strain (float): Strain value specified by the voigt index, vix.
            """
            super().__init__(np.zeros([3, 3]))
            self.set_by_voigt(vix, strain)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor)class DeformTensor(_Tensor3x3):
        """Deformation tensor object."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.DeformTensor.html#matlantis_features.features.elasticity.tensors.DeformTensor.__init__)    def __init__(self, strain_matrix: StrainTensor) -> None:
            """Initialize an instance.
    
            Args:
                strain_matrix (StrainTensor):
                  Strain tensor object converted to the deformation tensor.
            """
            # Create deformation gradient tensor from Green-Lagrange strain tensor,
            # based on the definition of Green-Lagrange strain tensor.
            # See: https://github.com/materialsproject/pymatgen/blob/v2022.0.4/pymatgen/analysis/elasticity/strain.py#L244  # NOQA
            mat = scipy.linalg.sqrtm(strain_matrix.to_3x3_ndarray() * 2 + np.eye(3))
            super().__init__(mat)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor)class StressTensor(_TensorVoigt):
        """Stress tensor object."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StressTensor.html#matlantis_features.features.elasticity.tensors.StressTensor.__init__)    def __init__(self, stress: np.ndarray) -> None:
            """Initialize an instance.
    
            Args:
                stress (np.ndarray): Stress tensor value in Voigt's notation.
            """
            super().__init__(stress)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList)class StrainTensorList:
        """Tensor list object for the strain tensors."""
    
        diagonal_strains: List[float]
        off_diagonal_strains: List[float]
        data: List[List[StrainTensor]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList.__init__)    def __init__(self, data: List[List[StrainTensor]]) -> None:
            """Initialize an instance.
    
            Args:
                data (list[list[StrainTensor]]): 2D-list of the strain tensors.
            """
            if len(data) != 6:
                ValueError(f"invalid len(data): {len(data)}")
            self.data = data
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList.from_strains)    @classmethod
        def from_strains(cls, diagonal_strains: List[float], off_diagonal_strains: List[float]) -> "StrainTensorList":
            """Create StrainTensorList from the list of diagonal and off-diagonal strain values.
    
            Args:
                diagonal_strains (list[float]): list of diagonal strain values.
                off_diagonal_strains (list[float]): list of diagonal strain values.
            Returns:
                StrainTensorList : StrainTensorList object.
            """
            data = [[StrainTensor(vix, s) for s in select(is_voigt_diag(vix), diagonal_strains, off_diagonal_strains)] for vix in range(0, 6)]
            obj = cls(data)
            obj.diagonal_strains = diagonal_strains
            obj.off_diagonal_strains = off_diagonal_strains
            return obj
    
    
    
        @property
        def n_total_samples(self) -> int:
            """Number of total sample points.
    
            Returns:
                int : Number of total sample points.
            """
            return sum([len(i) for i in self.data])
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StrainTensorList.html#matlantis_features.features.elasticity.tensors.StrainTensorList.get_samples)    def get_samples(self, vix: int) -> List[float]:
            """Get the strain sample list for the fitting calculation.
    
            Args:
                vix (int): Voigt index corresponding to the data-point axis.
            Returns:
                list[float] : List of stress values for the fitting calculation.
            """
            return [i.get_by_voigt(vix) for i in self.data[vix]]
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.DeformTensorList.html#matlantis_features.features.elasticity.tensors.DeformTensorList)class DeformTensorList:
        """Tensor list object for the deformation tensors."""
    
        data: List[List[DeformTensor]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.DeformTensorList.html#matlantis_features.features.elasticity.tensors.DeformTensorList.__init__)    def __init__(self, data: List[List[DeformTensor]]) -> None:
            """Initialize an instance.
    
            Args:
                data (list[list[DeformTensor]]): 2D-list of the deformation tensors.
            """
            if len(data) != 6:
                ValueError(f"invalid len(data): {len(data)}")
            self.data = data
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.DeformTensorList.html#matlantis_features.features.elasticity.tensors.DeformTensorList.from_strain_tensors)    @classmethod
        def from_strain_tensors(cls, strain_list: StrainTensorList) -> "DeformTensorList":
            """Create deform tensor list from the strain tensor list.
    
            Args:
                strain_list (StrainTensorList):
                  Strain tensor list object converted to the deform tensor list.
            Returns:
                DeformTensorList : Deform tensor list object.
            """
            data = map2d(strain_list.data, DeformTensor)
            return cls(data)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StressTensorList.html#matlantis_features.features.elasticity.tensors.StressTensorList)class StressTensorList:
        """Tensor list object for the stress."""
    
        data: List[List[StressTensor]]
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StressTensorList.html#matlantis_features.features.elasticity.tensors.StressTensorList.__init__)    def __init__(self, data: List[List[StressTensor]]) -> None:
            """Initialize an instance.
    
            Args:
                data (list[list[StressTensor]]): 2D-list of the stress tensors.
            """
            if len(data) != 6:
                ValueError(f"invalid len(data): {len(data)}")
            self.data = data
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.tensors.StressTensorList.html#matlantis_features.features.elasticity.tensors.StressTensorList.get_samples)    def get_samples(self, vi: int, vj: int) -> List[float]:
            """Get the stress sample list for the fitting calculation.
    
            Args:
                vi (int): Voigt index corresponding to the data-point axis.
                vj (int): Voigt index corresponding to the stress tensor axis.
            Returns:
                list[float] : List of stress values for the fitting calculation.
            """
            return [i.get_by_voigt(vj) for i in self.data[vi]]
    
    
    
