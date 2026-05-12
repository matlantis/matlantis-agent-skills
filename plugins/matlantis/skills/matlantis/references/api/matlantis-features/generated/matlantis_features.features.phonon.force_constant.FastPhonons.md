# matlantis_features.features.phonon.force_constant.FastPhonons#

_class _matlantis_features.features.phonon.force_constant.FastPhonons(_calculator : Calculator_, _supercell : Tuple[int, int, int]_, _delta : float = 0.1_, _method : str = 'frederiksen'_, _direction : str = 'central'_, _symmetrize : int = 3_, _acoustic : bool = True_, _use_symmetry : bool = False_, _cutoff : Optional[float] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_)[[source]](../_modules/matlantis_features/features/phonon/force_constant.html#FastPhonons)#
    

Bases: `object`

Class for calculating force constant using the finite displacement method.

Methods

`__init__`(calculator, supercell[, delta, ...]) | Initialize an instance.  
---|---  
`run`(atoms) | Runs the force constant calculation.  
  
__init__(_calculator : Calculator_, _supercell : Tuple[int, int, int]_, _delta : float = 0.1_, _method : str = 'frederiksen'_, _direction : str = 'central'_, _symmetrize : int = 3_, _acoustic : bool = True_, _use_symmetry : bool = False_, _cutoff : Optional[float] = None_, _show_progress_bar : bool = False_, _tqdm_options : Optional[Dict[str, Any]] = None_, _show_logger : bool = False_)[[source]](../_modules/matlantis_features/features/phonon/force_constant.html#FastPhonons.__init__)#
    

Initialize an instance.

Parameters
    

  * **calculator** (_Calculator_) – The calculator to estimate energy and force.

  * **supercell** (_tuple_ _[__int_ _,__int_ _,__int_ _]_) – The number of repetitions of input structure.

  * **delta** (_float_ _,__optional_) – The distance of finite displacment. Defaults to 0.1.

  * **method** (_str_ _,__optional_) – The method to calculate force constant. Defaults to ‘frederiksen’.

  * **direction** (_str_ _,__optional_) – How the displacement is found. Defaults to ‘central’.

  * **symmetrize** (_int_ _,__optional_) – The number of interactions to symmetrize the force constant. Defaults to 3.

  * **acoustic** (_bool_ _,__optional_) – Restores the acoustic sum rule of the force constant. Defaults to True.

  * **use_symmetry** (_bool_) – If True, the symmetry of input structure will be used to reduce the number of force estimations. Defaults to False.

  * **cutoff** (_float_ _or_ _None_ _,__optional_) – The corresponding elements in the force constant is set to 0 if the interatomic distance between atoms is larger than cutoff. Defaults to None.

  * **show_progress_bar** (_bool_ _,__optional_) – Show progress bar. Defaults to False.

  * **tqdm_options** (_dict_ _[__str_ _,__Any_ _] or_ _None_ _,__optional_) – Options for tqdm.

  * **show_logger** (_bool_ _,__optional_) – Show log information. Defaults to False.




run(_atoms : Union[[MatlantisAtoms](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms"), Atoms]_) → ndarray[[source]](../_modules/matlantis_features/features/phonon/force_constant.html#FastPhonons.run)#
    

Runs the force constant calculation.

Parameters
    

**atoms** ([_MatlantisAtoms_](../matlantis_features.atoms.html#matlantis_features.atoms.MatlantisAtoms "matlantis_features.atoms.MatlantisAtoms") _or_ _ASEAtoms_) – The input structure.

Returns
    

The force constant

Return type
    

np.ndarray

[ __ previous matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult ](matlantis_features.features.phonon.force_constant.ForceConstantFeatureResult.html "previous page") [ next matlantis_features.features.phonon.mode __](matlantis_features.features.phonon.mode.html "next page")
