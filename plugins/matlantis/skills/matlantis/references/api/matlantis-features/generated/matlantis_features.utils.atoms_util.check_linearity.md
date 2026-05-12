# matlantis_features.utils.atoms_util.check_linearity#

matlantis_features.utils.atoms_util.check_linearity(_atoms : Atoms_, _eigen_tolerance : float = 0.0001_) → [MoleculeLinearity](matlantis_features.utils.atoms_util.MoleculeLinearity.html#matlantis_features.utils.atoms_util.MoleculeLinearity "matlantis_features.utils.atoms_util.MoleculeLinearity")[[source]](../_modules/matlantis_features/utils/atoms_util.html#check_linearity)#
    

Check the linearity of molecule. The linearity can be linear, nonlinear and monatomic.

Parameters
    

  * **atoms** (_Atoms_) – The input structure.

  * **eigen_tolerance** (_float_ _,__optional_) – Tolerance of eigen value of inertia tensor when counting the number of zero elements. Defaults to 0.0001.



Returns
    

The linearity of the molecule.

Return type
    

[MoleculeLinearity](matlantis_features.utils.atoms_util.MoleculeLinearity.html#matlantis_features.utils.atoms_util.MoleculeLinearity "matlantis_features.utils.atoms_util.MoleculeLinearity")

[ __ previous matlantis_features.utils.atoms_util.MoleculeLinearity ](matlantis_features.utils.atoms_util.MoleculeLinearity.html "previous page") [ next matlantis_features.utils.atoms_util.check_symmetry_number __](matlantis_features.utils.atoms_util.check_symmetry_number.html "next page")
