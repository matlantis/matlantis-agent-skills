# Source code for matlantis_features.features.elasticity.elastic_tensor
    
    
    import itertools
    import logging
    from functools import partial
    from typing import Any, Callable, Dict, List, Optional, Tuple, Union
    
    import numpy as np
    import scipy.stats
    from ase import Atoms, units
    from tqdm.auto import tqdm
    from tqdm.contrib.logging import logging_redirect_tqdm
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    from matlantis_features.features.base import FeatureBase
    from matlantis_features.features.common import ASEOptFeature
    from matlantis_features.utils import FeatureCost, get_calculator
    from matlantis_features.utils.calculators.estimator_fn import EstimatorFnType
    
    from .tensors import (
        DeformTensor,
        DeformTensorList,
        StrainTensorList,
        StressTensor,
        StressTensorList,
    )
    from .utils import (
        equal_elements,
        get_symmetry_lattice,
        inverse_elements,
        is_voigt_diag,
        map2d,
        zeros_elements,
    )
    
    logger = logging.getLogger(__name__)
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms.html#matlantis_features.features.elasticity.elastic_tensor.make_deformed_atoms)def make_deformed_atoms(atoms: Atoms, dmat: DeformTensor) -> Atoms:
        """Create deformed Atoms object by deformation tensor.
    
        Args:
            atoms (Atoms): Input atoms.
            dmat (DeformTensor): Deformation tensor.
        Returns:
            Atoms : Resulting Atoms object.
        """
        cell = np.array(atoms.get_cell())
        deformed_cell = cell.dot(dmat.to_3x3_ndarray())
        atoms_deform = atoms.copy()
        atoms_deform.set_cell(deformed_cell, scale_atoms=True)
        return atoms_deform
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.stress_from_deform.html#matlantis_features.features.elasticity.elastic_tensor.stress_from_deform)def stress_from_deform(
        dmat: DeformTensor,
        atoms: Atoms,
        opt_fn: Callable[[Atoms], Atoms],
        show_logger: bool,
        pbar: Optional[tqdm],
        estimator_fn: Optional[EstimatorFnType] = None,
    ) -> StressTensor:
        """Calculate the stress tensor from the target atoms and deformation tensor, \
        using specified optimization method.
    
        Args:
            dmat (DeformTensor): Deformation tensor applied to the system.
            atoms (Atoms): Atoms representing the target system.
            opt_fn (Atoms -> Atoms): Optimization function
                to optimize the system after the deformation.
            show_logger (bool): Show log information.
            pbar (tqdm or None): Progress bar.
            estimator_fn (EstimatorFnType or None, optional):
                A factory method to create a custom estimator.
                Please refer :ref:`custom_estimator` for detail.
        Returns:
            StressTensor : Resulting stress tensor.
        """
        with logging_redirect_tqdm():
            if show_logger:
                logger.info("Calculating stress tensor")
                logger.info("Deform " + np.array2string(dmat.to_voigt_ndarray(), formatter={"float_kind": lambda x: "%.6f" % x}))
        atoms_deform = make_deformed_atoms(atoms, dmat)
        atoms_deform_opt = opt_fn(atoms_deform)
        with get_calculator(estimator_fn=estimator_fn) as calculator:
            atoms_deform_opt.calc = calculator
            stress = atoms_deform_opt.get_stress(voigt=True)
    
        with logging_redirect_tqdm():
            if show_logger:
                logger.info("Stress " + np.array2string(stress, formatter={"float_kind": lambda x: "%.2e" % x}))
            if pbar is not None:
                pbar.update(1)
    
        return StressTensor(stress)
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor)class ElasticTensor:
        """Elastic tensor object calculated from the stress-strain relation."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.__init__)    def __init__(self, stress_data: np.ndarray, strain_data: np.ndarray) -> None:
            """Initialize an instance.
    
            Args:
                stress_data (np.ndarray): Stress tensor for the calculation.
                strain_data (np.ndarray): Strain tensor for the calculation.
            """
            self.stress_data = stress_data
            self.strain_data = strain_data
    
            self.elastic_tensor = np.zeros([6, 6])
            self.intercepts = np.zeros([6, 6])
            self.rvalues = np.zeros([6, 6])
            self.pvalues = np.zeros([6, 6])
            self.std_errs = np.zeros([6, 6])
    
            self._calc_elastic_tensor()
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.from_stress_strain)    @classmethod
        def from_stress_strain(cls, stress: StressTensorList, strain: StrainTensorList) -> "ElasticTensor":
            """Create elastic tensor object from stress and strain tensor objects.
    
            Args:
                stress (StressTensorList): Stress tensor list object.
                strain (StrainTensorList): Strain tensor list object.
            Returns:
                ElasticTensor : Resulting elastic tensor object.
            """
            stress_data = np.empty(shape=(6, 6), dtype=object)
            strain_data = np.empty(shape=(6, 6), dtype=object)
    
            for i in range(6):
                strain_i = strain.get_samples(i)
                if not is_voigt_diag(i):
                    strain_i = [x * 2 for x in strain_i]
                for j in range(6):
                    stress_j = stress.get_samples(i, j)
                    logger.debug(f"strain: {strain_i}")
                    logger.debug(f"stress: {stress_j}")
                    if len(strain_i) != len(stress_j):
                        raise ValueError(f"inconsistent strain/strain sample sizes {len(strain_i)} != {len(stress_j)}")
                    stress_data[i, j] = stress_j
                    strain_data[i, j] = strain_i
            return cls(stress_data, strain_data)
    
    
    
        def _calc_elastic_tensor(self) -> None:
            """Calculate elastic tensor from stress and strain, assuming their linear relationship.
    
            This code is based on:
                Scientific Data (2), 150009 (2015) DOI: 10.1038/sdata.2015.9
            """
            for i, j in itertools.product(range(6), range(6)):
                logger.debug(f"solve {i}, {j}")
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(self.strain_data[i, j], self.stress_data[i, j])
                self.elastic_tensor[i, j] = slope
                self.intercepts[i, j] = intercept
                self.rvalues[i, j] = r_value
                self.pvalues[i, j] = p_value
                self.std_errs[i, j] = std_err
    
            self.elastic_tensor *= 1.0 / units.GPa
            self.elastic_tensor = 0.5 * (self.elastic_tensor + self.elastic_tensor.transpose())
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensor.fitting_result_linear)    def fitting_result_linear(self, i: int, j: int, threshold: float = 0.5) -> bool:
            """Check the linear relation between stress and strain values.
    
            Args:
                i (int): i-index for the tensor.
                j (int): j-index for the tensor.
                threshold (float, optional): Linearlity criterion threshold. Defaults to 0.5.
            Returns:
                bool : Returns True if the relation satisfy the linearlity criterion threshold.
            """
            ec = self.elastic_tensor[i, j]
            if abs(ec) < 1e-7:
                logger.warning(f"{i} {j}: ec is too small")
                return False
            det = abs(self.std_errs[i, j] / ec)
            if det > threshold:
                logger.debug(f"{i} {j} ec: {ec} err: {self.std_errs[i, j]} det: {det}")
                return False
            return True
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature)class ElasticTensorFeature(FeatureBase):
        """The matlantis-feature for calculating the elastic tensor for the specified atoms, \
        by applying the specified strain.
    
        References:
            For more information, please refer to `Matlantis Guidebook <https://redirect.matlantis.com/documents/detail/documents/matlantis-guidebook/en/matlantis-features/elastic_tensor.html>`_.
        """
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.__init__)    def __init__(
            self,
            diagonal_strains: List[float],
            off_diagonal_strains: List[float],
            optimizer: ASEOptFeature,
            check_linear_elems: Optional[List[Tuple[int, int]]] = None,
            zero_tol: float = 1e-3,
            show_progress_bar: bool = False,
            tqdm_options: Optional[Dict[str, Any]] = None,
            show_logger: bool = False,
            check_symmetry: bool = False,
            rtol: float = 0.01,
            atol: float = 1.0,
            estimator_fn: Optional[EstimatorFnType] = None,
        ):
            """Initialize an instance.
    
            Args:
                diagonal_strains (list[float]): Diagonal elements of the strain tensor.
                off_diagonal_strains (list[float]): Off-diagonal elements of the strain tensor.
                optimizer (ASEOptFeature): Optimization method to optimize the system after applying the strain.
                check_linear_elems (list[tuple[int, int]] or None, optional):
                  Index list of the elastic tensor to check the linearlity.
                  Default index list ((0, 0),(1, 1),(2, 2),(0, 1),(0, 2),(1, 2),(3, 3),(4, 4),(5, 5)) is used,
                  if None is specified.
                  Defaults to None.
                zero_tol (float, optional):
                  Skip the linearity check if the elastic constant is smaller than this threshold.
                  Defaults to 0.001.
                show_progress_bar (bool, optional): Show progress bar. Defaults to False.
                tqdm_options (dict[str, Any] or None, optional): Options for tqdm.
                show_logger (bool, optional): Show log information. Defaults to False.
                check_symmetry (bool, optional):
                  Check the calculated elastic tensor is in accord with the type of lattice symmetry.
                  For example, C11 = C22 and C14 = 0 in the cubic lattice. Please note that
                  the two elements a and b are considered as equal if their difference is within a tolerance
                  that us defined by rtol and atol, i.e. absolute(a - b) <= (atol + rtol * absolute(b)).
                  If check_symmetry is enabled, the resultant elastic tensor corresponds to the rotated atomic coordinates
                  after applying the symmetry operations, which might differ from the input atomic coordinates.
                  Please call get_symmetry_lattice(atoms) if the rotated atomic coordinates are needed.
                  In other words, if you do not want the atomic coordinates to be rotated from the original,
                  please set `check_symmetry` to False. Default to False.
                rtol (float, optional):
                  The relative tolerance of difference. The elastic tensor elements are regarded as the same
                  if the relative difference is smaller than this value. Default to 0.01.
                atol (float, optional):
                  The absolute tolerance of difference. The elastic tensor elements are regarded as the same
                  if the absolute difference is smaller than this value. Default to 1.0
                estimator_fn (EstimatorFnType or None, optional):
                    A factory method to create a custom estimator.
                    Please refer :ref:`custom_estimator` for detail.
            """
            super().__init__()
            self.check_estimator_fn(estimator_fn=estimator_fn)
            with self.init_scope():
                self.strain_list = StrainTensorList.from_strains(diagonal_strains, off_diagonal_strains)
                self.n_samples = self.strain_list.n_total_samples
                logger.debug(f"N: {self.n_samples}")
                logger.debug(f"self.strain_mats: {self.strain_list.data}")
    
                self.opt = optimizer.repeat(self.n_samples)
                self.show_progress_bar = show_progress_bar
                self.tqdm_options = tqdm_options
                self.show_logger = show_logger
    
                self.zero_tol = zero_tol
                if check_linear_elems is None:
                    self.check_linear_elems = [
                        (0, 0),
                        (1, 1),
                        (2, 2),
                        (0, 1),
                        (0, 2),
                        (1, 2),
                        (3, 3),
                        (4, 4),
                        (5, 5),
                    ]
                    # The check linear elems depends on the lattice symmetry
                    # If the triclinic lattice, all elements of elastic tensor are non-zero.
                else:
                    self.check_linear_elems = check_linear_elems
    
                self.check_symmetry = check_symmetry
                self.rtol = rtol
                self.atol = atol
    
                if self.opt.filter:
                    logging.warning("The cell optimization is enable in the 'optimizer', and this will lead to unrealistic result.")
    
                self.estimator_fn = estimator_fn
    
    
    
        def __call__(self, matlantis_atoms: Union[ASEAtoms, MatlantisAtoms]) -> ElasticTensor:
            """Calculate the elastic tensor.
    
            Args:
                matlantis_atoms (ASEAtoms or MatlantisAtoms):
                  Atoms representing the system to calculate the elastic tensor.
            Returns:
                ElasticTensor : Resulting elastic tensor object.
            """
            atoms: ASEAtoms
            if isinstance(matlantis_atoms, MatlantisAtoms):
                atoms = matlantis_atoms.ase_atoms
            else:
                atoms = matlantis_atoms
    
            if self.check_symmetry:
                atoms = get_symmetry_lattice(atoms)
    
            pbar = None
            if self.show_progress_bar:
                ntotal = self.n_samples
                if self.tqdm_options is None:
                    pbar = tqdm(total=ntotal)
                else:
                    pbar = tqdm(total=ntotal, **self.tqdm_options)
    
            deform_list = DeformTensorList.from_strain_tensors(self.strain_list)
            logger.debug(f"deform_list: {deform_list.data}")
    
            def opt_fn(atoms_deform: Atoms) -> Atoms:
                deform_opt = self.opt(atoms_deform)
                return deform_opt.atoms.ase_atoms
    
            stress_list = StressTensorList(
                map2d(
                    deform_list.data,
                    partial(
                        stress_from_deform,
                        atoms=atoms,
                        opt_fn=opt_fn,
                        show_logger=self.show_logger,
                        pbar=pbar,
                        estimator_fn=self.estimator_fn,
                    ),
                )
            )
            if pbar is not None:
                pbar.close()
    
            elastic_tensor = ElasticTensor.from_stress_strain(stress_list, self.strain_list)
            if self.check_symmetry:
                self.check_elastic_tensor_symmetry(atoms, elastic_tensor)
    
            return elastic_tensor
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.check_elastic_tensor_symmetry)    def check_elastic_tensor_symmetry(self, atoms: ASEAtoms, elastic_tensor: ElasticTensor) -> None:
            """Check correctness of elastic tensor according to the lattice symmetry.
    
            1. Certain elements of elastic tensor should be zeros, except the triclinic lattice.
            2. Some elements of elastic tensor should be identical, e.g. C11 = C22 = C33 in cubic lattice.
            3. Some elements of elastic tensor should be opposite number, e.g. C16 = -C26 in tetragonal lattice.
            4. Certain relation exists between the elements, e.g. C66 = 1/2 (C11-C12) in hexagonal lattice.
    
            Args:
                atoms (ASEAtoms):
                  Atoms representing the system to calculate the elastic tensor.
                elastic_tensor (ElasticTensor):
                  Elastic tensor object containing the stress and strain data.
            """
            rtol = self.rtol
            atol = self.atol
            cell = atoms.get_cell()
            lattice = cell.get_bravais_lattice()
            tensor = elastic_tensor.elastic_tensor
            lattice_system = lattice.lattice_system
            for i, j in zeros_elements[lattice_system]:
                if not np.isclose(tensor[i, j], 0.0, rtol=rtol, atol=atol):
                    logging.warning(
                        f"Warning: Elastic tensor of {lattice_system} lattice. C{i + 1}{j + 1} ({tensor[i, j]:2.3f}) is not close to 0 ."
                    )
            for (i1, j1), (i2, j2) in equal_elements[lattice_system]:
                if not np.isclose(tensor[i1, j1], tensor[i2, j2], rtol=rtol, atol=atol):
                    logging.warning(
                        f"Warning: Elastic tensor of {lattice_system} lattice. "
                        f"C{i1 + 1}{j1 + 1} ({tensor[i1, j1]:2.3f}) is not close to "
                        f"C{i2 + 1}{j2 + 1} ({tensor[i2, j2]:2.3f}). "
                    )
            for (i1, j1), (i2, j2) in inverse_elements[lattice_system]:
                if not np.isclose(tensor[i1, j1], -tensor[i2, j2], rtol=rtol, atol=atol):
                    logging.warning(
                        f"Warning: Elastic tensor of {lattice_system} lattice. "
                        f"C{i1 + 1}{j1 + 1} ({tensor[i1, j1]:2.3f}) is not close to the opposite of "
                        f"C{i2 + 1}{j2 + 1} ({tensor[i2, j2]:2.3f}). "
                    )
            if lattice_system in ["hexagonal", "Trigonal"]:
                if not np.isclose(tensor[5, 5], (tensor[0, 0] - tensor[0, 1]) / 2.0, rtol=rtol, atol=atol):
                    logging.warning(
                        f"Warning: Elastic tensor of {lattice_system} lattice. "
                        f"C66 ({tensor[5, 5]}) is not close to "
                        f"the (C11-C12)/2.0 ({(tensor[0, 0] - tensor[0, 1]) / 2.0})"
                    )
            for i in range(6):
                for j in range(i, 6):
                    if (np.abs(elastic_tensor.elastic_tensor[i, j]) > self.atol) and (not elastic_tensor.fitting_result_linear(i, j)):
                        msg = f"Nonlinear stress-strain relation in elastic tensor C{i + 1}{j + 1}"
                        logging.warning(msg)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.check_nonlinear)    def check_nonlinear(self, elastic_tensor: ElasticTensor, check_linear_elems: List[Tuple[int, int]]) -> None:
            """Check non-linearity of stress-strain relation.
    
            If non linearity is detected, warning will be reported to the logger
            and self.report() will be called.
    
            Args:
                elastic_tensor (ElasticTensor):
                  Elastic tensor object containing the stress and strain data.
                check_linear_elems (list[tuple[int, int]]):
                  Index list of the tensor elements to check the linearity.
            """
            for i, j in check_linear_elems:
                if not elastic_tensor.fitting_result_linear(i, j):
                    msg = f"Nonlinear stress-strain relation in elastic tensor C{i + 1}{j + 1}"
                    logging.warning(msg)
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.html#matlantis_features.features.elasticity.elastic_tensor.ElasticTensorFeature.check_zero_elems)    def check_zero_elems(
            self,
            elastic_tensor: ElasticTensor,
            check_linear_elems: List[Tuple[int, int]],
            zero_tol: float,
        ) -> None:
            """Check if the elastic tensor elements are zero.
    
            If non zero emenets are detected, warning will be reported to the logger
            and self.report() will be called.
    
            Args:
                elastic_tensor (ElasticTensor):
                  Elastic tensor object to be checked.
                check_linear_elems (list[tuple[int, int]]):
                  Index list of the tensor elements to check the linearity.
                  The elements specified in this list will not be checked,
                  since these elements are expected to have non-zero elastic constants.
                zero_tol (float): Tolerance threshold for the zero check.
            """
            for i, j in itertools.product(range(6), range(6)):
                if (i, j) not in check_linear_elems and (j, i) not in check_linear_elems:
                    ec = elastic_tensor.elastic_tensor[i, j]
                    if abs(ec) > zero_tol:
                        msg = f"Non-zero elastic tensor C{i}{j}={ec}"
                        logging.warning(msg)
    
    
    
        def _internal_cost(self, n_atoms: int) -> FeatureCost:
            """A function used for keeping track of cost information.
    
            Args:
                n_atoms (int): Number of atoms used for the calculation.
            Returns:
                FeatureCost : Cost information for the feature.
            """
            return FeatureCost(n_pfp_run=0, time_additional=1.0)
    
    
    
