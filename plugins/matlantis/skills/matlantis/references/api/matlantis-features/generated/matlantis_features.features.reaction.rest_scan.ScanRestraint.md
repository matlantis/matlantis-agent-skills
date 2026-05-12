# matlantis_features.features.reaction.rest_scan.ScanRestraint#

_class _matlantis_features.features.reaction.rest_scan.ScanRestraint(_direction : int_, _destination : float_, _exceed : float_, _fmax : float_)[[source]](../_modules/restscan/restraints.html#ScanRestraint)#
    

Bases: `object`

Super class for scan restraints.

Parameters
    

  * **direction** – 1 for scan towords large value. -1 for scan for towards small value.

  * **destination** – The scan destination.

  * **exceed** – In this class, scan force is turned off when it reaches the destination. This is a value indicating how far the scan can exceed the destination. This value is defined to smooth out the force field that Optimizer feels.

  * **fmax** – Rsstraint scan force.




Methods

`__init__`(direction, destination, exceed, fmax) |   
---|---  
`adjust_forces`(atoms, forces) |   
`adjust_momenta`(atoms, momenta) |   
`adjust_positions`(atoms, positions) |   
`get_coefficient`(atoms) |   
`get_colvar`(atoms) |   
`get_colvar_forces`(atoms) |   
`get_removed_dof`(atoms) |   
  
__init__(_direction : int_, _destination : float_, _exceed : float_, _fmax : float_) → None[[source]](../_modules/restscan/restraints.html#ScanRestraint.__init__)#
    

adjust_forces(_atoms : Atoms_, _forces : ndarray_) → None[[source]](../_modules/restscan/restraints.html#ScanRestraint.adjust_forces)#
    

adjust_momenta(_atoms : Atoms_, _momenta : ndarray_) → None[[source]](../_modules/restscan/restraints.html#ScanRestraint.adjust_momenta)#
    

adjust_positions(_atoms : Atoms_, _positions : ndarray_) → None[[source]](../_modules/restscan/restraints.html#ScanRestraint.adjust_positions)#
    

get_coefficient(_atoms : Atoms_) → float[[source]](../_modules/restscan/restraints.html#ScanRestraint.get_coefficient)#
    

get_colvar(_atoms : Atoms_) → float[[source]](../_modules/restscan/restraints.html#ScanRestraint.get_colvar)#
    

get_colvar_forces(_atoms : Atoms_) → ndarray[[source]](../_modules/restscan/restraints.html#ScanRestraint.get_colvar_forces)#
    

get_removed_dof(_atoms : Atoms_) → int[[source]](../_modules/restscan/restraints.html#ScanRestraint.get_removed_dof)#
    

[ __ previous matlantis_features.features.reaction.rest_scan.ScanDistance ](matlantis_features.features.reaction.rest_scan.ScanDistance.html "previous page") [ next matlantis_features.features.reaction.rest_scan.calculate_dihedral __](matlantis_features.features.reaction.rest_scan.calculate_dihedral.html "next page")
