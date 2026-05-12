# matlantis_features.features.phonon.utils.KptPath#

_class _matlantis_features.features.phonon.utils.KptPath(_cell : Cell_, _labels : Optional[List[str]] = None_, _special_kpts : Optional[ndarray] = None_, _n_kpts : Optional[List[int]] = None_, _total_n_kpts : int = 100_, _style : str = 'ase'_)[[source]](../_modules/matlantis_features/features/phonon/utils.html#KptPath)#
    

Bases: `object`

Class for k-point band path.

Methods

`__init__`(cell[, labels, special_kpts, ...]) | Creates the k-point band path.  
---|---  
  
__init__(_cell : Cell_, _labels : Optional[List[str]] = None_, _special_kpts : Optional[ndarray] = None_, _n_kpts : Optional[List[int]] = None_, _total_n_kpts : int = 100_, _style : str = 'ase'_)[[source]](../_modules/matlantis_features/features/phonon/utils.html#KptPath.__init__)#
    

Creates the k-point band path.

Parameters
    

  * **cell** (_Cell_) – The cell shape of input structure.

  * **labels** (_list_ _[__str_ _] or_ _None_ _,__optional_) – The labels of special points that form the band path. Please use ‘,’ ‘|’ or ‘/’ to represent a discontinuous jump. If None is provided, the band path will be automatically generated according to the type of bravais lattice. Defaults to None.

  * **special_kpts** (_np.ndarray_ _or_ _None_ _,__optional_) – The coordinates of special points in reciprocal space. The ‘special_kpts’ should be a numpy.ndarray in the shape of (N, 3), where N is the same as the length of ‘labels’. If ‘None’ is provided, the automatically generated special points will be used. Defaults to None.

  * **n_kpts** (_list_ _[__int_ _] or_ _None_ _,__optional_) – The number of interpolated points in each segment of the path. If None is provided, the values will be estimated from ‘total_n_kpts’. Defaults to None.

  * **total_n_kpts** (_int_ _,__optional_) – The number of interpolated points along the whole band path. This parameter will be ignored if ‘n_kpts’ is specified. Defaults to 100.

  * **style** (_str_ _,__optional_) – If ‘labels’ and ‘special_kpts’ are all not provided, the default k-point path will be used. The parameter ‘style’ can control which kind of default path to use. Currently, ‘ase’, which is ASE style, and ‘phonon_db’, which follows the Phonon database (<http://phonondb.mtl.kyoto-u.ac.jp/>), are supported. Defaults to “ase”.




[ __ previous matlantis_features.features.phonon.utils ](matlantis_features.features.phonon.utils.html "previous page") [ next matlantis_features.features.phonon.utils.PhononFrequency __](matlantis_features.features.phonon.utils.PhononFrequency.html "next page")
