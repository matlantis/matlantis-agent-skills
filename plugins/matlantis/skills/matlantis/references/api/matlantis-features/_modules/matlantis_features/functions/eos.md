# Source code for matlantis_features.functions.eos
    
    
    import warnings
    from typing import Callable, List, Tuple, Union
    
    import numpy as np
    from scipy.optimize import curve_fit, fsolve
    
    
    # TODO add other types of EOS: Birch Murnaghan, Vinet ...
    class EOSFitError(Exception):
        pass
    
    
    class EOSVolumeError(Exception):
        pass
    
    
    def _initial_guess(V: Union[np.ndarray, List[float]], E: Union[np.ndarray, List[float]]) -> Tuple[float, float, float, float]:
        coeff: np.ndarray = np.polyfit(V, E, 2)
        V0: float = np.abs(-coeff[1] / 2 / coeff[0])
        E0: float = coeff[2] + coeff[1] * V0 + coeff[0] * V0**2
        B0: float = np.abs(2 * coeff[0] * V0)
        return E0, B0, 4.0, V0
    
    
    def _birch_murnaghan_func(V: Union[np.ndarray, float], E0: float, B0: float, BP: float, V0: float) -> Union[np.ndarray, float]:
        eta: Union[np.ndarray, float] = (V0 / V) ** (2 / 3.0)
        E = E0 + 9 * B0 * V0 / 16 * (eta - 1) ** 2 * (6 + BP * (eta - 1) - 4 * eta)
        return E
    
    
    def _birch_murnaghan_func_deriv(V: Union[np.ndarray, float], E0: float, B0: float, BP: float, V0: float) -> Union[np.ndarray, float]:
        eta: Union[np.ndarray, float] = (V0 / V) ** (2 / 3.0)
        dE_dV = -3 * V0 * B0 / 8 * eta * (eta - 1) * (3 * (eta - 1) * BP - 12 * eta + 16) / V
        return dE_dV
    
    
    
    
    [[docs]](../../../generated/matlantis_features.functions.eos.BirchMurnaghanEOS.html#matlantis_features.functions.eos.BirchMurnaghanEOS)class BirchMurnaghanEOS:
    
    
    [[docs]](../../../generated/matlantis_features.functions.eos.BirchMurnaghanEOS.html#matlantis_features.functions.eos.BirchMurnaghanEOS.__init__)    def __init__(self, E0: float, B0: float, BP: float, V0: float) -> None:
            """The Birch–Murnaghan equation of state.
    
            Args:
                E0 (float): The energy at the equilibrium volume.
                B0 (float): The bulk modulus.
                BP (float): the derivative of the bulk modulus with respect to pressure.
                V0 (float): The equilibirum volume.
            """
            self.E0 = E0
            self.B0 = B0
            self.BP = BP
            self.V0 = V0
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.functions.eos.BirchMurnaghanEOS.html#matlantis_features.functions.eos.BirchMurnaghanEOS.__call__)    def __call__(self, V: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
            """Estimate the energy at the given volume from the Birch–Murnaghan equation of state.
    
            Args:
                V (np.ndarray): Volume.
            Returns:
                np.ndarray : Energy.
            """
            return _birch_murnaghan_func(V, self.E0, self.B0, self.BP, self.V0)
    
    
    
    
    
    [[docs]](../../../generated/matlantis_features.functions.eos.BirchMurnaghanEOS.html#matlantis_features.functions.eos.BirchMurnaghanEOS.from_fit)    @classmethod
        def from_fit(
            cls,
            V: Union[np.ndarray, List[float]],
            E: Union[np.ndarray, List[float]],
            max_nfev: int = 100000,
        ) -> "BirchMurnaghanEOS":
            """Get the Birch–Murnaghan equation of state from fitting a set of volumes and energies.
    
            Args:
                V (Union[np.ndarray, List[float]]): Volumes.
                E (Union[np.ndarray, List[float]]): Energies.
                max_nfev (int): The maximum number of function calls to fit the curve.
            Returns:
                BirchMurnaghanEOS : The fitted Birch–Murnaghan equation of state.
            """
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                try:
                    E0, B0, BP, V0 = _initial_guess(V, E)
                    popt, pcov = curve_fit(
                        _birch_murnaghan_func,
                        V,
                        E,
                        p0=[E0, B0, BP, V0],
                        max_nfev=max_nfev,
                        bounds=([-np.inf, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
                    )
                except (ValueError, RuntimeError, TypeError, Warning) as e:
                    print(e)
                    raise EOSFitError
                else:
                    return cls(E0=popt[0], B0=popt[1], BP=popt[2], V0=popt[3])
    
    
    
        @property
        def bulk_modulus(self) -> float:
            """The bulk modulus of the material estimated from the Birch–Murnaghan equation of state.
    
            Returns:
                float : The bulk modulus.
            """
            return self.B0
    
        @property
        def equilibrium_volume(self) -> float:
            """The eqilibrium volume of the material estimated from the Birch–Murnaghan equation of state.
    
            Returns:
                float : The eqilibrium volume.
            """
            return self.V0
    
    
    
    [[docs]](../../../generated/matlantis_features.functions.eos.BirchMurnaghanEOS.html#matlantis_features.functions.eos.BirchMurnaghanEOS.get_volume_under_pressure)    def get_volume_under_pressure(self, pressure: float) -> float:
            """Get the volume under a given pressure.
    
            Args:
                pressure (float): The external pressure.
            Returns:
                float : The volume.
            """
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
    
                    def func(
                        c: float,
                    ) -> Callable[[Union[np.ndarray, float]], Union[np.ndarray, float]]:
                        def func_(V: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
                            V_p = _birch_murnaghan_func_deriv(V, self.E0, self.B0, self.BP, self.V0) + c
                            return V_p
    
                        return func_
    
                    V: float = fsolve(func(pressure), self.V0)[0]
                except ValueError:
                    raise EOSVolumeError
                else:
                    return V
    
    
    
