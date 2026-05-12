# matlantis_features.functions.eos.BirchMurnaghanEOS#

_class _matlantis_features.functions.eos.BirchMurnaghanEOS(_E0 : float_, _B0 : float_, _BP : float_, _V0 : float_)[[source]](../_modules/matlantis_features/functions/eos.html#BirchMurnaghanEOS)#
    

Bases: `object`

Methods

`__init__`(E0, B0, BP, V0) | The Birch–Murnaghan equation of state.  
---|---  
`__call__`(V) | Estimate the energy at the given volume from the Birch–Murnaghan equation of state.  
`from_fit`(V, E[, max_nfev]) | Get the Birch–Murnaghan equation of state from fitting a set of volumes and energies.  
`get_volume_under_pressure`(pressure) | Get the volume under a given pressure.  
  
__init__(_E0 : float_, _B0 : float_, _BP : float_, _V0 : float_) → None[[source]](../_modules/matlantis_features/functions/eos.html#BirchMurnaghanEOS.__init__)#
    

The Birch–Murnaghan equation of state.

Parameters
    

  * **E0** (_float_) – The energy at the equilibrium volume.

  * **B0** (_float_) – The bulk modulus.

  * **BP** (_float_) – the derivative of the bulk modulus with respect to pressure.

  * **V0** (_float_) – The equilibirum volume.




__call__(_V : Union[ndarray, float]_) → Union[ndarray, float][[source]](../_modules/matlantis_features/functions/eos.html#BirchMurnaghanEOS.__call__)#
    

Estimate the energy at the given volume from the Birch–Murnaghan equation of state.

Parameters
    

**V** (_np.ndarray_) – Volume.

Returns
    

Energy.

Return type
    

np.ndarray

_classmethod _from_fit(_V : Union[ndarray, List[float]]_, _E : Union[ndarray, List[float]]_, _max_nfev : int = 100000_) → BirchMurnaghanEOS[[source]](../_modules/matlantis_features/functions/eos.html#BirchMurnaghanEOS.from_fit)#
    

Get the Birch–Murnaghan equation of state from fitting a set of volumes and energies.

Parameters
    

  * **V** (_Union_ _[__np.ndarray_ _,__List_ _[__float_ _]__]_) – Volumes.

  * **E** (_Union_ _[__np.ndarray_ _,__List_ _[__float_ _]__]_) – Energies.

  * **max_nfev** (_int_) – The maximum number of function calls to fit the curve.



Returns
    

The fitted Birch–Murnaghan equation of state.

Return type
    

BirchMurnaghanEOS

get_volume_under_pressure(_pressure : float_) → float[[source]](../_modules/matlantis_features/functions/eos.html#BirchMurnaghanEOS.get_volume_under_pressure)#
    

Get the volume under a given pressure.

Parameters
    

**pressure** (_float_) – The external pressure.

Returns
    

The volume.

Return type
    

float

[ __ previous matlantis_features.functions.eos ](matlantis_features.functions.eos.html "previous page") [ next Tips for matlantis-features __](../tips/index.html "next page")
