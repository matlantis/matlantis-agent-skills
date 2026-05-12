# matlantis_features.utils.atoms_util.convert_atoms_to_upper#

matlantis_features.utils.atoms_util.convert_atoms_to_upper(_atoms : Atoms_) → Atoms[[source]](../_modules/matlantis_features/utils/atoms_util.html#convert_atoms_to_upper)#
    

Rotate atoms ase NPT module requires atoms to have following upper-triangular form for cell: [[X X X] [0 X X] [0 0 X]] this method rotates atoms to make its cell upper triangular.

Parameters
    

**atoms** (_Atoms_) – input atoms.

Returns
    

output, rotated atoms.

Return type
    

rotated_atoms (Atoms)

[ __ previous matlantis_features.utils.atoms_util.check_symmetry_number ](matlantis_features.utils.atoms_util.check_symmetry_number.html "previous page") [ next matlantis_features.utils.atoms_util.get_inertia_tensor __](matlantis_features.utils.atoms_util.get_inertia_tensor.html "next page")
