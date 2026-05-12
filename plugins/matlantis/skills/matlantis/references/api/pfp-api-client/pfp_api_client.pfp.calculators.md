# Calculator interface#

_class _pfp_api_client.pfp.calculators.ase_calculator.ASECalculator(_estimator : [Estimator](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator "pfp_api_client.pfp.estimator.Estimator")_)[[source]](_modules/pfp_api_client/pfp/calculators/ase_calculator.html#ASECalculator)#
    

Bases: `Calculator`

ASECalculator is a derivation of ase.calculators.calculator. Usually this will be used from ase.Atoms object. It can be safe to have multiple instances of ASECalculator running concurrently in different threads, as long as they do not share the same Estimator instance. It is _not_ safe to use a single instance of ASECalculator concurrently and trying to do so will likely raise a ConcurrentUseDetected exception.

Parameters
    

**estimator** ([_Estimator_](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator "pfp_api_client.pfp.estimator.Estimator")) – PFP Estimator

atoms _: Optional[Atoms]_#
    

calculate(_atoms : Optional[Atoms] = None_, _properties : Optional[Sequence[str]] = None_, _system_changes : Sequence[str] = ['positions', 'numbers', 'cell', 'pbc', 'initial_charges', 'initial_magmoms']_) → None[[source]](_modules/pfp_api_client/pfp/calculators/ase_calculator.html#ASECalculator.calculate)#
    

Do the calculation. This function will be called from various ASE objects such as ase.optimizer or ase.md classes.

Parameters
    

  * **atoms** – ase.Atoms object.

  * **properties** – Sequence of what needs to be calculated. Current supported values are ‘energy’, ‘forces’, ‘stress’, and ‘charges’.

  * **system_changes** – Sequence of what has changed since last calculation.




current_messages _: Set[[MessageEnum](pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum "pfp_api_client.pfp.utils.messages.MessageEnum")]_#
    

estimator _: [Estimator](pfp_api_client.pfp.estimator.html#pfp_api_client.pfp.estimator.Estimator "pfp_api_client.pfp.estimator.Estimator")_#
    

get_descriptors(_atoms : Atoms_, _indices : Optional[ndarray] = None_) → Dict[str, ndarray][[source]](_modules/pfp_api_client/pfp/calculators/ase_calculator.html#ASECalculator.get_descriptors)#
    

Calculate and return scalar descriptors for a given set of atoms.

Parameters
    

  * **atoms** (_ase.Atoms_) – The ASE Atoms object for which descriptors are to be calculated.

  * **indices** (_Optional_ _[__np.ndarray_ _]__,__optional_) – A 1-dimensional array of indices (between 0 and n_atoms-1) specifying the atomic indices for descriptors to be obtained.



Returns
    

**descriptors** – Dictionary containing the calculated scalar descriptors under the key scalar. The descriptors have a shape of [n_atoms, 256].

Return type
    

Dict[str, np.ndarray]

set_default_properties(_default_properties : Sequence[str]_) → None[[source]](_modules/pfp_api_client/pfp/calculators/ase_calculator.html#ASECalculator.set_default_properties)#
    

Specify the properties that are always calculated.

Parameters
    

**default_properties** (_Sequence_ _[__str_ _]_) – Properties to always calculate. Must be a subset of `implemented_properties`.
