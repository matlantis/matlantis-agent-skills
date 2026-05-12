  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.convex...



# mtcsp.convex_hull_search package#

## Module contents#

### Experiment configuration#

_class _mtcsp.convex_hull_search.BaseSearchConfig(_elements_ , _experiment_name_ , _pfp_version =None_, _pfp_method_type =None_, _pfp_compatibility =Compatibility.PFP2023COMPATIBILITY_)#
    

Dataclass for common search configurations for convex hull search and derivative structure search.

Parameters:
    

  * **elements** (_Elements_)

  * **experiment_name** (_str_)

  * **pfp_version** (_InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_method_type** (_InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_compatibility** (_InitVar_ _[__Compatibility_ _|__str_ _]_)




elements _: `Annotated`[`tuple`[`str`, `...`]]_#
    

List of elements to search for.

experiment_name _: `str`_#
    

Unique name for the experiment.

pfp_compatibility _: `InitVar`_ _ = 'pfp2023compatibility'_#
    

Scheme to choose and combine PFP’s calc_mode across entire composition space.

There are two options:

  * Compatibility.PFP2023COMPATIBILITY: (Default) Mixes EstimatorCalcMode.PBE and EstimatorCalcMode.PBE_U depending on compositions, which corresponds to the GGA/GGA+U and anion corrections in Materials Project.

  * Compatibility.R2SCAN: Use EstimatorCalcMode.R2SCAN for all compositions.

  * Compatibility.R2SCAN_PLUS_D3: (Provisional) Use EstimatorCalcMode.R2SCAN_PLUS_D3 for all compositions.

  * Compatibility.PBE: Use EstimatorCalcMode.PBE for all compositions.

  * Compatibility.PBE_PLUS_D3: (Provisional) Use EstimatorCalcMode.PBE_PLUS_D3 for all compositions.




pfp_method_type _: `InitVar`_ _ = None_#
    

Method type of the PFP to use. If None, the default method type will be used.

If you would like to use MNCore, please set ‘mncore’ here.

pfp_version _: `InitVar`_ _ = None_#
    

Version of the PFP to use. If None, the latest version on Matlantis will be used.

If you would like to use a consistent version, please provide it here.

_class _mtcsp.convex_hull_search.SearchConfig(_elements_ , _experiment_name_ , _pfp_version =None_, _pfp_method_type =None_, _pfp_compatibility =Compatibility.PFP2023COMPATIBILITY_, _population_size =128_, _max_atoms =32_, _add_mp_crystals_to_initial_population =False_, _composition_range =None_, _composition_endpoints =None_, _oxidation_state_ranges =None_, _seed =None_)#
    

Bases: `BaseSearchConfig`

Dataclass for configuration of single search.

Parameters:
    

  * **elements** (_Annotated_ _[__tuple_ _[__str_ _,__...__]__,__'Tuple_ _of_ _chemical symbols'__]_)

  * **experiment_name** (_str_)

  * **pfp_version** (_dataclasses.InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_method_type** (_dataclasses.InitVar_ _[__str_ _|__None_ _]_)

  * **pfp_compatibility** (_dataclasses.InitVar_ _[__mtcsp.localopt._pfp.Compatibility_ _|__str_ _]_)

  * **population_size** (_int_)

  * **max_atoms** (_int_)

  * **add_mp_crystals_to_initial_population** (_bool_)

  * **composition_range** (_dict_ _[__str_ _,__tuple_ _[__Fraction_ _,__Fraction_ _]__]__|__None_)

  * **composition_endpoints** (_list_ _[__str_ _|_[_Element_](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Element "\(in pymatgen v2023.7.17\)") _|_[_Species_](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Species "\(in pymatgen v2023.7.17\)") _|_[_DummySpecies_](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.DummySpecies "\(in pymatgen v2023.7.17\)") _|__dict_ _|_[_Composition_](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition "\(in pymatgen v2023.7.17\)") _]__|__None_)

  * **oxidation_state_ranges** (_dict_ _[__str_ _,__tuple_ _[__int_ _,__int_ _]__]__|__None_)

  * **seed** (_int_ _|__None_)




add_mp_crystals_to_initial_population _: `bool`_ _ = False_#
    

Whether to add MP crystals to the initial population.

composition_endpoints _: `list`[`str` | [`Element`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Element "\(in pymatgen v2023.7.17\)") | [`Species`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Species "\(in pymatgen v2023.7.17\)") | [`DummySpecies`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.DummySpecies "\(in pymatgen v2023.7.17\)") | `dict` | [`Composition`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition "\(in pymatgen v2023.7.17\)")] | `None`_ _ = None_#
    

If provided, the search is limited to the convex hull defined by these endpoints.

Each endpoint is a composition (e.g., “SrO”, “SrO2”, “TiO2”) representing a vertex of the simplex (line, triangle, etc.) explored during the search. All compositions within the convex region spanned by these endpoints are considered.

composition_range _: `dict`[`str`, `tuple`[`Fraction`, `Fraction`]] | `None`_ _ = None_#
    

If provided, the search is limited to this composition.

Each key-value pair specifies an element and tuple of two fractions.Fraction that represents the permissible composition range (minimum and maximum). Every element should be specified in the composition_range even if the range is the full range of 0 to 1 (e.g., {“O”: (Fraction(0), Fraction(1))}). The composition range for each element is applied independently.

For example, composition_range={“Ti”: (Fraction(“1/3”), Fraction(“1/2”)), “O”: (Fraction(0), Fraction(1))} limits the normalized composition to Ti_{x}O_{y} with constraints 1/3 <= x <= 1/2, 0 <= y <=1, and x + y = 1. The additional constraint x + y = 1 is assumed by default and does not need to be explicitly specified.

Please avoid using float to ensure the endpoints are included in the permissible range.

max_atoms _: `int`_ _ = 32_#
    

Maximum number of atoms in structures.

oxidation_state_ranges _: `dict`[`str`, `tuple`[`int`, `int`]] | `None`_ _ = None_#
    

If provided, the search is limited to compositions satisfying these oxidation state ranges.

Each key-value pair specifies an element and a tuple of two integers that represent the minimum and maximum values of the permissible oxidation state, respectively.

For example, oxidation_state_ranges={“Ti”: (2, 4), “O”: (-2, -2)} limits the oxidation states of Ti to be between +2 and +4 (inclusive) and O to be -2.

For single element and alloy-like compositions, the oxidation state ranges are ignored.

population_size _: `int`_ _ = 128_#
    

Population size for the search.

seed _: `int` | `None`_ _ = None_#
    

If provided, the random seed for the search will be set to this value.

validate()#
    

Validate the search configuration. :meta private:

Return type:
    

`None`

_class _mtcsp.convex_hull_search.SystemConfig(_db_file_ , _atoms_store_dir_ , _atoms_store_type ='filesystem'_, _log_file =None_)#
    

Configuration for the system.

Parameters:
    

  * **db_file** (_Path_)

  * **atoms_store_dir** (_Path_)

  * **atoms_store_type** (_Literal_ _[__'filesystem'__,__'zipfiles'__]_)

  * **log_file** (_Path_ _|__None_)




atoms_store_dir _: `Path`_#
    

Directory to store Atoms files.

atoms_store_type _: `Literal`[`'filesystem'`, `'zipfiles'`]__ = 'filesystem'_#
    

Type of Atoms store.

db_file _: `Path`_#
    

Path for database. If storage ends with “.journal”, JournalStorage is used, otherwise RDBStorage is used.

_classmethod _from_experiment_name(_experiment_name_ , _*_ , _atoms_store_type ='filesystem'_)#
    

Return type:
    

`Self`

Parameters:
    

  * **experiment_name** (_str_)

  * **atoms_store_type** (_Literal_ _[__'filesystem'__,__'zipfiles'__]_)




log_file _: `Path` | `None`_ _ = None_#
    

Path for log file. If None, directly print to stderr.

### Experiment initialization and execution#

mtcsp.convex_hull_search.initialize_hull_search(_search_config_ , _system_config_)#
    

Initialize and return an experiment for convex hull search.

Parameters:
    

  * **search_config** (`SearchConfig`) – Configuration for the convex hull search.

  * **system_config** (`SystemConfig`) – System configuration.



Return type:
    

[`Experiment`](mtcsp.experiment.html#mtcsp.experiment.Experiment "mtcsp.experiment._experiment.Experiment")

Returns:
    

Initialized experiment for convex hull search.

mtcsp.convex_hull_search.perform_hull_search(_*_ , _search_config_ , _system_config_ , _n_crystal_structures_ , _parallelism =6_, _n_jobs_per_sampler_process =4_)#
    

Perform convex hull search.

Parameters:
    

  * **search_config** (`SearchConfig`) – Configuration for the convex hull search.

  * **system_config** (`SystemConfig`) – System configuration.

  * **n_crystal_structures** (`int`) – Number of crystal structures to search. The search will stop when the number of successfully completed trials reaches this number.

  * **parallelism** (`int`) – Number of parallel processes to use, by default 6. (default: `6`)

  * **n_jobs_per_sampler_process** (_int_)



Return type:
    

`None`

[ previous mtcsp.atoms package ](mtcsp.atoms.html "previous page") [ next mtcsp.crystal_structure package ](mtcsp.crystal_structure.html "next page")

On this page 

  * Module contents
    * Experiment configuration
      * `BaseSearchConfig`
        * `BaseSearchConfig.elements`
        * `BaseSearchConfig.experiment_name`
        * `BaseSearchConfig.pfp_compatibility`
        * `BaseSearchConfig.pfp_method_type`
        * `BaseSearchConfig.pfp_version`
      * `SearchConfig`
        * `SearchConfig.add_mp_crystals_to_initial_population`
        * `SearchConfig.composition_endpoints`
        * `SearchConfig.composition_range`
        * `SearchConfig.max_atoms`
        * `SearchConfig.oxidation_state_ranges`
        * `SearchConfig.population_size`
        * `SearchConfig.seed`
        * `SearchConfig.validate()`
      * `SystemConfig`
        * `SystemConfig.atoms_store_dir`
        * `SystemConfig.atoms_store_type`
        * `SystemConfig.db_file`
        * `SystemConfig.from_experiment_name()`
        * `SystemConfig.log_file`
    * Experiment initialization and execution
      * `initialize_hull_search()`
      * `perform_hull_search()`


  *[*]: Keyword-only parameters separator (PEP 3102)
