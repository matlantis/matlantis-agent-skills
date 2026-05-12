  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.deriva...



# mtcsp.derivative_structure package#

## Module contents#

### Experiment configuration#

_class _mtcsp.derivative_structure.DerivativeStructureSearchConfig(_elements_ , _experiment_name_ , _pfp_version =None_, _pfp_method_type =None_, _pfp_compatibility =Compatibility.PFP2023COMPATIBILITY_)#
    

Bases: [`BaseSearchConfig`](mtcsp.convex_hull_search.html#mtcsp.convex_hull_search.BaseSearchConfig "mtcsp.convex_hull_search._config.BaseSearchConfig")

Dataclass for configuration of single derivative structure search.

Parameters:
    

  * **elements** (_Annotated_ _[__tuple_ _[__str_ _,__...__]__,__'Tuple_ _of_ _chemical symbols'__]_)

  * **experiment_name** (_str_)

  * **pfp_version** (_dataclasses.InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_method_type** (_dataclasses.InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_compatibility** (_dataclasses.InitVar_ _[__mtcsp.localopt._pfp.Compatibility_ _|__str_ _]_)




_classmethod _from_feasible_number_mapping(_feasible_number_mapping_ , _experiment_name_ , _pfp_version =None_, _pfp_method_type =None_, _pfp_compatibility =Compatibility.PFP2023COMPATIBILITY_)#
    

Create a `DerivativeStructureSearchConfig` from feasible_number_mapping.

Parameters:
    

  * **feasible_number_mapping** (`dict`[`int`, `list`[`int`]]) – Mapping from atomic number to feasible atomic numbers in the derivative structures.

  * **experiment_name** (`str`) – Name of the experiment.

  * **pfp_version** (`str` | `None`) – Version of the PFP model to use. If None, the default version will be used. (default: `None`)

  * **pfp_method_type** (`str` | `None`) – Method type of the PFP model to use. If None, the default method type will be used. (default: `None`)

  * **pfp_compatibility** (`Compatibility` | `str`) – Compatibility scheme to use for the PFP model. Default is Compatibility.get_default(). (default: `<Compatibility.PFP2023COMPATIBILITY: 'pfp2023compatibility'>`)



Return type:
    

`Self`

### Derivative structure generation#

mtcsp.derivative_structure.append_derivative_structures(_experiment_ , _base_structure_ , _max_index_)#
    

Append derivative structures to the experiment.

Symmetrically distinct derivative structures up to the given max_index are enumerated and prepended to the experiment.

Parameters:
    

  * **experiment** (`DerivativeStructureExperiment`) – Derivative structure experiment to append structures to.

  * **base_structure** (`BaseStructure`) – Base structure to enumerate derivative structures from.

  * **max_index** (`int`) – Maximum supercell size (index) to enumerate derivative structures. Use `list_num_derivative_structures()` to determine appropriate value.



Return type:
    

`None`

mtcsp.derivative_structure.append_derivative_structures_with_random_sampling(_experiment_ , _base_structure_ , _supercell_ , _n_crystal_structures_ , _seed =0_)#
    

Append derivative structures to the experiment using random sampling.

Symmetrically equivalent derivative structures may be repeatedly sampled.

Parameters:
    

  * **experiment** (`DerivativeStructureExperiment`) – Derivative structure experiment to append structures to.

  * **base_structure** (`BaseStructure`) – Base structure to enumerate derivative structures from.

  * **supercell** (`list`[`list`[`int`]]) – Supercell matrix to generate derivative structures. An obtained derivative structure’s j-th axis is given by sum_i supercell[i][j] * base_structure.basis[i].

  * **n_crystal_structures** (`int`) – Number of derivative structures to sample and append.

  * **seed** (`int`) – Random seed for sampling. Default is 0. (default: `0`)



Return type:
    

`None`

### Experiment initialization and execution#

mtcsp.derivative_structure.initialize_derivative_structure_search(_search_config_ , _system_config_)#
    

Initialize and return derivative structure search.

Parameters:
    

  * **search_config** (`DerivativeStructureSearchConfig`) – Configuration for the derivative structure search.

  * **system_config** ([`SystemConfig`](mtcsp.convex_hull_search.html#mtcsp.convex_hull_search.SystemConfig "mtcsp.convex_hull_search._config.SystemConfig")) – Configuration for the system.



Return type:
    

`DerivativeStructureExperiment`

Returns:
    

Initialized derivative structure experiment.

mtcsp.derivative_structure.perform_derivative_structure_search(_*_ , _search_config_ , _system_config_ , _parallelism =6_)#
    

Perform derivative structure search.

Parameters:
    

  * **search_config** (`DerivativeStructureSearchConfig`) – Configuration for the derivative structure search.

  * **system_config** ([`SystemConfig`](mtcsp.convex_hull_search.html#mtcsp.convex_hull_search.SystemConfig "mtcsp.convex_hull_search._config.SystemConfig")) – Configuration for the system.

  * **parallelism** (`int`) – Number of parallel processes to use. Default is 6. (default: `6`)



Return type:
    

`None`

### Utility functions#

mtcsp.derivative_structure.prepare_base_structure_from_atoms(_atoms_ , _feasible_number_mapping_ , _composition =None_)#
    

Prepare a BaseStructure from an ASE Atoms object.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – ASE Atoms object representing the base structure.

  * **feasible_number_mapping** (`dict`[`int`, `list`[`int`]]) – Mapping from atomic numbers to feasible atomic numbers. Let x be an atomic number in atoms.numbers, then a proceeding enumeration will create a structure by replacing atomic species of x with one of feasible_number_mapping[x].

  * **composition** (`dict`[`int`, `int`] | `None`) – 

Composition of the **enumerated** structures. Specify the integer ratio of atoms for each feasible atomic number. The integer ratios are not required to be mutually coprime.

For example, to enumerate LiCo0.5Ni0.5O2 from the base structure LiCoO2, use the following:
        
        composition={
            ase.data.atomic_numbers['Li']: 2,
            ase.data.atomic_numbers['Co']: 1,
            ase.data.atomic_numbers['Ni']: 1,
            ase.data.atomic_numbers['O']: 4,
        } (default: ``None``)
        



Return type:
    

`BaseStructure`

mtcsp.derivative_structure.list_num_derivative_structures(_base_structure_ , _*_ , _max_total_num_structures =30000_)#
    

Return a list of the number of symmetrically distinct derivative structures for increasing supercell sizes (index).

Parameters:
    

  * **base_structure** (`BaseStructure`) – Base structure to enumerate derivative structures from.

  * **max_total_num_structures** (`int`) – Maximum total number of derivative structures to enumerate. The enumeration continues until the total number of derivative structures reaches this value. (default: `30000`)



Return type:
    

`list`[`dict`]

Returns:
    

A list of dictionaries containing the following keys:

  * ”index”: Supercell size

  * ”num_structures”: Number of symmetrically distinct derivative structures for the given index

  * ”total_num_structures”: Cumulative total number of derivative structures up to the given index




mtcsp.derivative_structure.VACANCY_NUMBER#
    

alias of 0

[ previous mtcsp.crystal_structure package ](mtcsp.crystal_structure.html "previous page") [ next mtcsp.experiment package ](mtcsp.experiment.html "next page")

On this page 

  * Module contents
    * Experiment configuration
      * `DerivativeStructureSearchConfig`
        * `DerivativeStructureSearchConfig.from_feasible_number_mapping()`
    * Derivative structure generation
      * `append_derivative_structures()`
      * `append_derivative_structures_with_random_sampling()`
    * Experiment initialization and execution
      * `initialize_derivative_structure_search()`
      * `perform_derivative_structure_search()`
    * Utility functions
      * `prepare_base_structure_from_atoms()`
      * `list_num_derivative_structures()`
      * `VACANCY_NUMBER`


  *[*]: Keyword-only parameters separator (PEP 3102)
