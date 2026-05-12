# matlantis_features.features.reaction.rest_scan.ScanDistance#

_class _matlantis_features.features.reaction.rest_scan.ScanDistance(_indices : Tuple[int, int]_, _direction : int_, _destination : float_, _exceed : float = 0.1_, _fmax : float = 0.05_)[[source]](../_modules/restscan/restraints.html#ScanDistance)#
    

Bases: [`ScanRestraint`](matlantis_features.features.reaction.rest_scan.ScanRestraint.html#matlantis_features.features.reaction.rest_scan.ScanRestraint "restscan.restraints.ScanRestraint")

Scan for distance.

Parameters
    

  * **direction** – 1 for scan towords large distance. -1 for scan for towards small distance.

  * **destination** – The scan finish when the distance reach this value.

  * **exceed** – In this class, scan force is turned off when it reaches the destination. This is a value indicating how far the scan can exceed the destination. This value is defined to smooth out the force field that Optimizer feels.

  * **fmax** – Rsstraint scan force.




Methods

`__init__`(indices, direction, destination[, ...]) |   
---|---  
`adjust_forces`(atoms, forces) |   
`adjust_momenta`(atoms, momenta) |   
`adjust_positions`(atoms, positions) |   
`get_coefficient`(atoms) |   
`get_colvar`(atoms) |   
`get_colvar_forces`(atoms) |   
`get_removed_dof`(atoms) |   
`index_shuffle`(atoms, ind) |   
`todict`() |   
  
__init__(_indices : Tuple[int, int]_, _direction : int_, _destination : float_, _exceed : float = 0.1_, _fmax : float = 0.05_) → None[[source]](../_modules/restscan/restraints.html#ScanDistance.__init__)#
    

adjust_forces(_atoms : Atoms_, _forces : ndarray_) → None#
    

adjust_momenta(_atoms : Atoms_, _momenta : ndarray_) → None#
    

adjust_positions(_atoms : Atoms_, _positions : ndarray_) → None#
    

get_coefficient(_atoms : Atoms_) → float#
    

get_colvar(_atoms : Atoms_) → float[[source]](../_modules/restscan/restraints.html#ScanDistance.get_colvar)#
    

get_colvar_forces(_atoms : Atoms_) → ndarray[[source]](../_modules/restscan/restraints.html#ScanDistance.get_colvar_forces)#
    

get_removed_dof(_atoms : Atoms_) → int#
    

index_shuffle(_atoms : Atoms_, _ind : ndarray_) → None[[source]](../_modules/restscan/restraints.html#ScanDistance.index_shuffle)#
    

todict() → Dict[str, Any][[source]](../_modules/restscan/restraints.html#ScanDistance.todict)#
    

[ __ previous matlantis_features.features.reaction.rest_scan.ScanDihedral ](matlantis_features.features.reaction.rest_scan.ScanDihedral.html "previous page") [ next matlantis_features.features.reaction.rest_scan.ScanRestraint __](matlantis_features.features.reaction.rest_scan.ScanRestraint.html "next page")
