# matlantis_features.features.phonon.utils.PhononFrequency#

_class _matlantis_features.features.phonon.utils.PhononFrequency(_force_constant : ndarray_, _supercell : ndarray_, _unit_cell_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_)[[source]](../_modules/matlantis_features/features/phonon/utils.html#PhononFrequency)#
    

Bases: `object`

The class to calculate the phonon frequencies and modes from the force constant.

Methods

`__init__`(force_constant, supercell, ...) | Initialize an instance.  
---|---  
`get_frequency`(kpts[, unit]) | Calculate the phonon frequency for each k-point.  
`get_frequency_mode`(kpts[, unit]) | Calculate both the phonon frequency and the phonon mode for each k-point.  
  
__init__(_force_constant : ndarray_, _supercell : ndarray_, _unit_cell_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_)[[source]](../_modules/matlantis_features/features/phonon/utils.html#PhononFrequency.__init__)#
    

Initialize an instance.

Parameters
    

  * **force_constant** (_np.ndarray_) – The force constant.

  * **supercell** (_np.ndarray_) – The supercell size used for the force constant calculation.

  * **unit_cell_atoms** ([_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_) – The input primitive cell structure.




get_frequency(_kpts : ndarray_, _unit : str = 'meV'_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/utils.html#PhononFrequency.get_frequency)#
    

Calculate the phonon frequency for each k-point.

Parameters
    

  * **kpts** (_ndarray_) – The phonon wavevectors (k-points), expected to be a 2D ndarray with shape (n,3).

  * **unit** (_str_ _,__optional_) – The unit of phonon frequency. Must be one of: ‘meV’, ‘eV’, ‘THz’ or ‘cm^-1’. Defaults to ‘meV’.



Returns
    

The phonon frequency.

Return type
    

ndarray

Raises
    

  * **ValueError** – If the kpts input is not a (n, 3) 2D ndarray.

  * **ValueError** – If the unit is not one of the allowed values.




get_frequency_mode(_kpts : ndarray_, _unit : str = 'meV'_) → Tuple[ndarray, ndarray][[source]](../_modules/matlantis_features/features/phonon/utils.html#PhononFrequency.get_frequency_mode)#
    

Calculate both the phonon frequency and the phonon mode for each k-point.

Parameters
    

  * **kpts** (_ndarray_) – The phonon wavevectors (k-points), expected to be a 2D ndarray with shape (n,3).

  * **unit** (_str_ _,__optional_) – The unit of phonon frequency. Must be one of: ‘meV’, ‘eV’, ‘THz’ or ‘cm^-1’. Defaults to ‘meV’.



Returns
    

The phonon frequency and the phonon mode.

Return type
    

tuple[ndarray, ndarray]

Raises
    

  * **ValueError** – If the kpts input is not a (n, 3) 2D ndarray.

  * **ValueError** – If the unit is not one of the allowed values.




[ __ previous matlantis_features.features.phonon.utils.KptPath ](matlantis_features.features.phonon.utils.KptPath.html "previous page") [ next matlantis_features.features.phonon.utils.get_k_mesh __](matlantis_features.features.phonon.utils.get_k_mesh.html "next page")
