# Source code for matlantis_features.features.phonon.utils
    
    
    import logging
    from typing import List, Optional, Tuple, Union
    
    import numpy as np
    import seekpath
    import spglib
    from ase import Atoms, units
    from ase.cell import Cell
    from ase.lattice import identify_lattice
    from numpy import ndarray
    
    from matlantis_features.atoms import ASEAtoms, MatlantisAtoms
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.PhononFrequency.html#matlantis_features.features.phonon.utils.PhononFrequency)class PhononFrequency(object):
        """The class to calculate the phonon frequencies and modes from the force constant."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.PhononFrequency.html#matlantis_features.features.phonon.utils.PhononFrequency.__init__)    def __init__(
            self,
            force_constant: np.ndarray,
            supercell: np.ndarray,
            unit_cell_atoms: Union[MatlantisAtoms, ASEAtoms],
        ):
            """Initialize an instance.
    
            Args:
                force_constant (np.ndarray): The force constant.
                supercell (np.ndarray): The supercell size used for the force constant
                    calculation.
                unit_cell_atoms (MatlantisAtoms or ASEAtoms): The input primitive cell structure.
            """
            if isinstance(unit_cell_atoms, MatlantisAtoms):
                unit_cell_ase_atoms: ASEAtoms = unit_cell_atoms.ase_atoms
            else:
                unit_cell_ase_atoms = unit_cell_atoms
    
            self.force_constant = force_constant
            self.supercell = supercell
            self.n_cells = np.prod(self.supercell)
            self.n_atoms = len(unit_cell_ase_atoms.get_atomic_numbers())
            self.numbers = unit_cell_ase_atoms.get_atomic_numbers()
            self.cell = unit_cell_ase_atoms.get_cell()
            self.mass = unit_cell_ase_atoms.get_masses()
    
    
    
        def _dynamic_matrix_q(self, q: ndarray) -> ndarray:
            """
            Calculate the dynamical matrix for a given q-point.
    
            This method calculates the dynamical matrix for phonon calculations
            at a given q-point, using the Fourier transform of the force constants
            and mass factors.
    
            Args:
                q (ndarray): The phonon wavevector (q-point), expected to be a 1D array with shape (3,).
    
            Returns:
                ndarray: The dynamical matrix corresponding to the given q-point.
    
            Raises:
                ValueError: If the q input is not a (3,) ndarray.
    
            References:
                - X. Gonze and C. Lee, Phys. Rev. B - Condens. Matter Mater. Phys. 55, 10355 (1997).
            """
            # check shape of q
            if not isinstance(q, ndarray) or q.shape != (3,):
                raise ValueError("Invalid input value for q-point. The q should be a (3,) ndarray.")
    
            # get the lattice vectors
            R = np.indices(self.supercell).reshape(3, -1)  # type: ignore
            for i in range(3):
                R[i][R[i] > self.supercell[i] // 2] -= self.supercell[i]
    
            # phase factor
            phase = np.exp(-2.0j * np.pi * np.dot(q, R))
            # print(phase)
    
            # Fourier transform
            C_q = np.sum(phase[:, np.newaxis, np.newaxis] * self.force_constant, axis=0)
    
            # Mass factor
            mass_factor = np.outer(np.repeat(self.mass**-0.5, 3), np.repeat(self.mass**-0.5, 3))
            D_q: ndarray = mass_factor * C_q
    
            # Force to be hermite matrix. This improves numerical stability of band structures.
            D_q = 0.5 * (D_q + D_q.conj().transpose())
    
            return D_q
    
        def _frequency_mode_q(self, q: ndarray, return_mode: bool = False) -> Union[Tuple[ndarray, ndarray], ndarray]:
            # get dynamic matrix
            dynamic_matrix = self._dynamic_matrix_q(q)
    
            # get eigenvalue (w) and eigenvector (u) of dynamic matrix
            if return_mode:
                w, u = np.linalg.eigh(dynamic_matrix, UPLO="U")
            else:
                w = np.linalg.eigvalsh(dynamic_matrix, UPLO="U")
    
            # vibraition frequency
            omega = np.empty_like(w)
            omega[w > 0] = np.sqrt(w[w > 0])
            omega[w < 0] = -np.sqrt(-w[w < 0])  # negtive frequency case
    
            # change the unit --> meV
            energy: ndarray = omega * units._hbar * 1e13 / np.sqrt(units._e * units._amu)  # units meV
    
            if return_mode:
                # vibition mode
                mode: ndarray = (np.repeat(self.mass**-0.5, 3)[:, np.newaxis] * u).T.reshape((-1, self.n_atoms, 3))
                return energy, mode
    
            else:
                return energy
    
        def _validate_kpts(self, kpts: ndarray) -> None:
            """Validate the inputs for kpts.
    
            Args:
                kpts (ndarray): The phonon wavevectors (k-points), expected to be a 2D ndarray with shape (n,3).
            """
            if not isinstance(kpts, ndarray) or kpts.ndim != 2 or kpts.shape[1] != 3:
                raise ValueError("Invalid input value for k-points. The k-points should be a (n, 3) 2D ndarray.")
    
        def _validate_unit(self, unit: str) -> None:
            """Validate the inputs for unit.
    
            Args:
                unit (str): The unit of phonon frequency. Must be one of: 'meV', 'eV', 'THz' or 'cm^-1'.
            """
            units = ["meV", "eV", "THz", "cm^-1"]
            if unit not in units:
                raise ValueError(f"Invalid unit value. Must be one of {units}.")
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.PhononFrequency.html#matlantis_features.features.phonon.utils.PhononFrequency.get_frequency)    def get_frequency(self, kpts: ndarray, unit: str = "meV") -> ndarray:
            """Calculate the phonon frequency for each k-point.
    
            Args:
                kpts (ndarray): The phonon wavevectors (k-points), expected to be a 2D ndarray with shape (n,3).
                unit (str, optional): The unit of phonon frequency. Must be one of: 'meV', 'eV', 'THz'
                    or 'cm^-1'. Defaults to 'meV'.
    
            Returns:
                ndarray : The phonon frequency.
    
            Raises:
                ValueError: If the kpts input is not a (n, 3) 2D ndarray.
                ValueError: If the unit is not one of the allowed values.
            """
            self._validate_kpts(kpts=kpts)
            self._validate_unit(unit=unit)
            self.n_kpts = len(kpts)
    
            energy_list = []
            for q in kpts:
                f = self._frequency_mode_q(q)
                energy_list.append(f)
            energy: ndarray = np.stack(energy_list)  # type: ignore
            if unit == "meV":
                pass
            elif unit == "eV":
                energy *= 1 / 1000.0
            elif unit == "THz":
                energy *= 0.24180
            else:
                energy *= 8.06558
    
            return energy
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.PhononFrequency.html#matlantis_features.features.phonon.utils.PhononFrequency.get_frequency_mode)    def get_frequency_mode(self, kpts: ndarray, unit: str = "meV") -> Tuple[ndarray, ndarray]:
            """Calculate both the phonon frequency and the phonon mode for each k-point.
    
            Args:
                kpts (ndarray): The phonon wavevectors (k-points), expected to be a 2D ndarray with shape (n,3).
                unit (str, optional): The unit of phonon frequency. Must be one of: 'meV', 'eV', 'THz'
                    or 'cm^-1'. Defaults to 'meV'.
            Returns:
                tuple[ndarray, ndarray] : The phonon frequency and the phonon mode.
    
            Raises:
                ValueError: If the kpts input is not a (n, 3) 2D ndarray.
                ValueError: If the unit is not one of the allowed values.
            """
            self._validate_kpts(kpts=kpts)
            self._validate_unit(unit=unit)
            self.n_kpts = len(kpts)
    
            energy_list = []
            mode_list = []
            for q in kpts:
                f, m = self._frequency_mode_q(q, return_mode=True)
                energy_list.append(f)
                mode_list.append(m.real)
            energy = np.stack(energy_list)
            mode = np.stack(mode_list)
    
            if unit == "meV":
                pass
            elif unit == "eV":
                energy *= 1 / 1000.0
            elif unit == "THz":
                energy *= 0.24180
            else:
                energy *= 8.06558
    
            return energy, mode
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.get_primitive_structure.html#matlantis_features.features.phonon.utils.get_primitive_structure)def get_primitive_structure(atoms: Atoms, primitive_matrix: Union[np.ndarray, str], prec: float = 1e-3) -> Atoms:
        """Transforms the input unit cell to the primitive cell.
    
        Args:
            atoms (Atoms): The input structure.
            primitive_matrix (Union[np.ndarray, str]): the primitive cell basis vectors or 'auto'.
                If a 3x3 array is provided, the structure will be wrapped into the primitive cell
                according to the user-defined axis. If 'auto', the primitive cell will be automatically
                defined according to its symmetry.
            prec (float, optional): Wrapper the input structure into the primitive cell with certain
                precision. Defaults to 0.001.
        Returns:
            Atoms : The primitive cell.
        """
        logger = logging.getLogger(__name__)
    
        if isinstance(primitive_matrix, ndarray):
            logger.info("Try to get primitive cell according to 'primitive_matrix'. \n")
            assert primitive_matrix.shape == (
                3,
                3,
            ), "Please provide the 'primitive_matrix' with correct shape."
            primitive_cell = get_primitive_structure_from_primitive_matrix(atoms, primitive_matrix, prec)[0]
        elif primitive_matrix == "auto":
            logger.info("Try to get primitive cell automatically. \n")
            primitive_cell = get_primitive_structure_auto(atoms, prec)
        else:
            raise ValueError(
                "Invalid input value for 'primitive_matrix'. The 'primitive_matrix' should be either a 3x3 array or a string 'auto'. "
            )
        logger.info(f"The primitive cell contains {len(primitive_cell)} atoms. \nThe lattice of primitive cell is {primitive_cell.cell}. \n")
        return primitive_cell
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.get_primitive_structure_auto.html#matlantis_features.features.phonon.utils.get_primitive_structure_auto)def get_primitive_structure_auto(atoms: Atoms, prec: float = 1e-3) -> Atoms:
        """Transforms the input cell to the primitive cell by automatically analyzing its symmetry.
    
        Args:
            atoms (Atoms): The input structure.
            prec (float, optional): Wrapper the input structure into the primitive cell with certain
                precision. Defaults to 0.001.
        Returns:
            Atoms : The primitive cell.
        """
        cell, poisitions, numbers = spglib.find_primitive((atoms.get_cell()[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers()))
        return Atoms(cell=cell, scaled_positions=poisitions, numbers=numbers, pbc=[True, True, True])
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix.html#matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix)def get_primitive_structure_from_primitive_matrix(
        atoms: Atoms, primitive_matrix: np.ndarray, prec: float = 1e-3
    ) -> Tuple[Atoms, np.ndarray]:
        """Transforms the input unit cell to the primitive cell according to the axis provided.
    
        Args:
            atoms (Atoms): The input structure.
            primitive_matrix (np.ndarray): the primitive cell basis vectors.
            prec (float, optional): Wrapper the input structure into the primitive cell with certain
                precision. Defaults to 0.001.
        Returns:
            tuple[Atoms, np.ndarray] : The primitive cell.
        """
        trimmed_lattice = np.dot(primitive_matrix.T, atoms.cell)
        n_wrapped_cells = 1.0 / np.linalg.det(primitive_matrix)
        positions_in_new_lattice = np.dot(atoms.get_scaled_positions(), np.linalg.inv(primitive_matrix).T)
        positions_in_new_lattice -= np.floor(positions_in_new_lattice)
    
        wrapped_positions_list: List[np.ndarray] = []
        map_s2p_list: List[int] = []
        map_p2s_list: List[int] = []
        for i, p in enumerate(positions_in_new_lattice):
            new_atom = True
            for j, wp in enumerate(wrapped_positions_list):
                dist = np.linalg.norm((p - wp).dot(trimmed_lattice))
                if dist < prec:
                    map_p2s_list.append(j)
                    new_atom = False
                    break
            if new_atom:
                wrapped_positions_list.append(p)
                map_p2s_list.append(len(wrapped_positions_list) - 1)
                map_s2p_list.append(i)
    
        wrapped_positions = np.array(wrapped_positions_list)
        map_s2p = np.array(map_s2p_list)
        map_p2s = np.array(map_p2s_list)
    
        n_wrapped_atoms = len(wrapped_positions)
        for i in range(n_wrapped_atoms):
            ind = map_p2s == i
            if np.sum(ind) != n_wrapped_cells:
                raise ValueError
            else:
                if len(np.unique(atoms.numbers[ind])) != 1:
                    raise ValueError
    
        primitive_atoms = Atoms(
            positions=atoms.positions[map_s2p],
            numbers=atoms.numbers[map_s2p],
            cell=trimmed_lattice,
            pbc=True,
        )
        return primitive_atoms, map_p2s
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.get_k_mesh.html#matlantis_features.features.phonon.utils.get_k_mesh)def get_k_mesh(size: np.ndarray, method: str = "mp", shift: Optional[np.ndarray] = None) -> np.ndarray:
        """Construct a uniform sampling of k-space of given size.
    
        Args:
            size (np.ndarray): The size of k-point mesh grid.
            method (str, optional): Create the k-point mesh grid with Monkhorst-Pack scheme ('mp") or
                Gamma centered scheme ('gamma'). Defaults to 'mp'.
            shift (np.ndarray or None, optional): An array of shape (3,) containing the shift of the
                mesh grid in the direction along the corresponding reciprocal axes. Defaults to None.
        Returns:
            np.ndarray : The k-point mesh.
        """
        # https://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
        assert len(size) == 3
        assert method in ["mp", "gamma"]
        if shift is None:
            shift = np.zeros(3)
        qpoints = []
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    if method == "mp":
                        qpoints.append(
                            [
                                (i + shift[0] + 0.5) / size[0],
                                (j + shift[1] + 0.5) / size[1],
                                (k + shift[2] + 0.5) / size[2],
                            ]
                        )
                    else:
                        qpoints.append(
                            [
                                (i + shift[0]) / size[0],
                                (j + shift[1]) / size[1],
                                (k + shift[2]) / size[2],
                            ]
                        )
        qpoints = np.array(qpoints)  # type: ignore
        return qpoints - np.rint(qpoints)  # type: ignore
    
    
    
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.KptPath.html#matlantis_features.features.phonon.utils.KptPath)class KptPath:
        """Class for k-point band path."""
    
    
    
    [[docs]](../../../../generated/matlantis_features.features.phonon.utils.KptPath.html#matlantis_features.features.phonon.utils.KptPath.__init__)    def __init__(
            self,
            cell: Cell,
            labels: Optional[List[str]] = None,
            special_kpts: Optional[np.ndarray] = None,
            n_kpts: Optional[List[int]] = None,
            total_n_kpts: int = 100,
            style: str = "ase",
        ):
            """Creates the k-point band path.
    
            Args:
                cell (Cell): The cell shape of input structure.
                labels (list[str] or None, optional): The labels of special points that form the band
                    path. Please use ',' '|' or '/' to represent a discontinuous jump. If None is
                    provided, the band path will be automatically generated according to the type of
                    bravais lattice. Defaults to None.
                special_kpts (np.ndarray or None, optional): The coordinates of special points in
                    reciprocal space. The 'special_kpts' should be a numpy.ndarray in the shape of
                    (N, 3), where N is the same as the length of 'labels'. If 'None' is provided,
                    the automatically generated special points will be used. Defaults to None.
                n_kpts (list[int] or None, optional): The number of interpolated points in each
                    segment of the path. If None is provided, the values will be estimated from
                    'total_n_kpts'. Defaults to None.
                total_n_kpts (int, optional): The number of interpolated points along the whole band
                    path. This parameter will be ignored if 'n_kpts' is specified. Defaults to 100.
                style (str, optional): If 'labels' and 'special_kpts' are all not provided,
                    the default k-point path will be used. The parameter 'style' can control which
                    kind of default path to use. Currently, 'ase', which is ASE style, and 'phonon_db',
                    which follows the Phonon database (http://phonondb.mtl.kyoto-u.ac.jp/), are
                    supported. Defaults to "ase".
            """
            assert style in ["ase", "phonon_db"]
    
            self.cell = cell
            if style == "phonon_db" and labels is None and special_kpts is None:
                skp = seekpath.get_explicit_k_path(
                    (cell.array, np.array([[0.0, 0.0, 0.0]]), np.array([1])),
                    with_time_reversal=True,
                    reference_distance=0.025,
                )
                self.kpts = skp["explicit_kpoints_rel"]
                self.coords = skp["explicit_kpoints_linearcoord"]
                self.labels: List[str] = []
                special_kpt_indices_: List[int] = []
                for i, lab in enumerate(skp["explicit_kpoints_labels"]):
                    if lab != "":
                        self.labels.append(lab)
                        special_kpt_indices_.append(i)
                self.special_kpt_indices = np.array(special_kpt_indices_)
                self.special_kpts = self.kpts[special_kpt_indices_]
            else:
                self.lattice = identify_lattice(self.cell)[0]
                if not labels:
                    labels = [i for i in self.lattice.bandpath().path]
    
                self.labels = []
                self.breakpoints: List[int] = []
                for label in labels:
                    if label in [",", "|", "/"]:
                        self.breakpoints.append(len(self.labels))
                    elif label in ["1", "2", "3"]:
                        self.labels[-1] += label
                    else:
                        self.labels.append(label)
    
                self.n_special_kpts = len(self.labels)
    
                if special_kpts is not None:
                    if not labels:
                        raise ValueError("Please provide labels")
    
                    if special_kpts.shape[0] != self.n_special_kpts:
                        raise ValueError("The number of labels is inconsistent with the number of special_kpts")
    
                    self.special_kpts = special_kpts
                else:
                    kpts_ = []
                    for point in self.labels:
                        if point not in self.lattice.special_point_names:
                            raise ValueError(
                                f"Illegal name of special k-point in {self.lattice.longname}. {point} is not included in ",
                                self.lattice.special_point_names,
                            )
                        else:
                            kpts_.append(self.lattice.get_special_points()[point])
                    self.special_kpts = np.array(kpts_)
    
                if n_kpts is None:
                    n_kpts = self._get_n_kpts_from_total(total_n_kpts)
    
                if len(n_kpts) < self.n_special_kpts - 1:
                    raise ValueError()
                else:
                    self.n_kpts: List[int] = n_kpts[: self.n_special_kpts - 1] + [0]
                    for i in self.breakpoints:
                        self.n_kpts[i - 1] = 0
    
                kpts, coords, special_kpt_indices = self._interp_path()
                self.kpts = kpts
                self.coords = coords
                self.special_kpt_indices = special_kpt_indices
    
    
    
        def _get_n_kpts_from_total(self, total_n_kpts: int) -> List[int]:
            """Estimates the number of interplotations in each segment according to the total number of \
                k-points.
    
            Args:
                total_n_kpts (int): The total number of k-points along the whole band path.
            Returns:
                list[int] : The number of interpolations in each segment.
            """
            dist_list = []
            for i in range(self.n_special_kpts - 1):
                if i + 1 in self.breakpoints:
                    dist_list.append(0.0)
                else:
                    dist_list.append(
                        np.linalg.norm(self.special_kpts[i + 1] - self.special_kpts[i])  # type: ignore
                    )
            dist = np.array(dist_list)
            n_kpts = list(np.rint(dist * total_n_kpts / dist.sum()).astype(int))
            return n_kpts
    
        def _interp_path(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
            """Linear interpolatation between special points.
    
            Returns:
                tuple[np.ndarray, np.ndarray, np.ndarray] :
                    The k-points, coordinates of k-points, and the locations of special points.
            """
            bp = [0] + self.breakpoints + [len(self.labels)]
            for i in range(len(bp) - 1):
                if i == 0:
                    kpts, coords, special_kpt_indices = self._interp_path_segment(bp[i], bp[i + 1])
                else:
                    kpts_, coords_, special_kpt_indices_ = self._interp_path_segment(bp[i], bp[i + 1])
                    coords_ += coords[-1]
                    special_kpt_indices_ += len(kpts)
                    kpts = np.vstack((kpts, kpts_))
                    coords = np.concatenate((coords, coords_))
                    special_kpt_indices = np.concatenate((special_kpt_indices, special_kpt_indices_))
            return kpts, coords, special_kpt_indices
    
        def _interp_path_segment(self, init: int, end: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
            """Linear interpolation in a segment.
    
            Args:
                init (int): The initial k-point.
                end (int): The end k-point.
            Returns:
                tuple[np.ndarray, np.ndarray, np.ndarray] :
                    The k-points, coordinates of k-points, and the locations of special points.
            """
            labels = self.labels[init:end]
            special_kpts = self.special_kpts[init:end]
            n_kpts = self.n_kpts[init:end]
    
            kpts_list = []
            for i in range(len(labels) - 1):
                x = np.linspace(0, 1, n_kpts[i], endpoint=False)[:, np.newaxis]
                kpts_ = special_kpts[i] * (1 - x) + special_kpts[i + 1] * x
                kpts_list.append(kpts_)
            kpts_list.append(special_kpts[-1])
            kpts = np.vstack(kpts_list)
    
            dist = np.linalg.norm(kpts[1:] - kpts[:-1], axis=1)
            coords = np.concatenate((np.array([0.0]), np.cumsum(dist)))
    
            special_kpt_indices = np.concatenate((np.array([0]), np.cumsum(n_kpts[:-1])))
            return kpts, coords, special_kpt_indices
    
    
    
