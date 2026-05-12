# matlantis_features.features.phonon.utils#

Classes

[`KptPath`](matlantis_features.features.phonon.utils.KptPath.html#matlantis_features.features.phonon.utils.KptPath "matlantis_features.features.phonon.utils.KptPath")(cell[, labels, special_kpts, ...]) | Class for k-point band path.  
---|---  
[`PhononFrequency`](matlantis_features.features.phonon.utils.PhononFrequency.html#matlantis_features.features.phonon.utils.PhononFrequency "matlantis_features.features.phonon.utils.PhononFrequency")(force_constant, supercell, ...) | The class to calculate the phonon frequencies and modes from the force constant.  
  
Functions

[`get_k_mesh`](matlantis_features.features.phonon.utils.get_k_mesh.html#matlantis_features.features.phonon.utils.get_k_mesh "matlantis_features.features.phonon.utils.get_k_mesh")(size[, method, shift]) | Construct a uniform sampling of k-space of given size.  
---|---  
[`get_primitive_structure`](matlantis_features.features.phonon.utils.get_primitive_structure.html#matlantis_features.features.phonon.utils.get_primitive_structure "matlantis_features.features.phonon.utils.get_primitive_structure")(atoms, primitive_matrix) | Transforms the input unit cell to the primitive cell.  
[`get_primitive_structure_auto`](matlantis_features.features.phonon.utils.get_primitive_structure_auto.html#matlantis_features.features.phonon.utils.get_primitive_structure_auto "matlantis_features.features.phonon.utils.get_primitive_structure_auto")(atoms[, prec]) | Transforms the input cell to the primitive cell by automatically analyzing its symmetry.  
[`get_primitive_structure_from_primitive_matrix`](matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix.html#matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix "matlantis_features.features.phonon.utils.get_primitive_structure_from_primitive_matrix")(...) | Transforms the input unit cell to the primitive cell according to the axis provided.  
  
[ __ previous matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram ](matlantis_features.features.phonon.qha_gibbs_free_energy.draw_phase_diagram.html "previous page") [ next matlantis_features.features.phonon.utils.KptPath __](matlantis_features.features.phonon.utils.KptPath.html "next page")
