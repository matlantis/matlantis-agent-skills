  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.crysta...



# mtcsp.crystal_structure package#

## Module contents#

_class _mtcsp.crystal_structure.CrystalStructure(_id_ , _generation=None_ , _unoptimized_crystal_structure=None_ , _chemical_formula=None_ , _operator_type=None_ , _origin=None_ , _parents=None_ , _steps=None_ , _fmax=None_ , _optimizer_type=None_ , _optimized_crystal_structure=None_ , _potential_energy=None_ , _formation_energy=None_ , _mixing_energy= <property object>_, _rejection_reason=None_ , _opt_stats=None_ , _minimum_neighbors_generation=None_ , _minimum_neighbors_value=None_ , _is_external=False_ , _cell_size=None_ , _metric_values=None_)#
    

Dataclass for a crystal structure under exploration or one that has already been explored.

Parameters:
    

  * **id** (_CrystalStructureID_)

  * **generation** (_int_ _|__None_)

  * **unoptimized_crystal_structure** (_str_ _|__None_)

  * **chemical_formula** (_str_ _|__None_)

  * **operator_type** (_str_ _|__None_)

  * **origin** (_str_ _|__None_)

  * **parents** (_tuple_ _[__CrystalStructureID_ _,__...__]__|__None_)

  * **steps** (_int_ _|__None_)

  * **fmax** (_float_ _|__None_)

  * **optimizer_type** (_OptimizerType_ _|__None_)

  * **optimized_crystal_structure** (_str_ _|__None_)

  * **potential_energy** (_float_ _|__None_)

  * **formation_energy** (_float_ _|__None_)

  * **mixing_energy** (_InitVar_ _[__float_ _|__None_ _]_)

  * **rejection_reason** (_str_ _|__None_)

  * **opt_stats** (_dict_ _[__str_ _,__float_ _|__None_ _]__|__None_)

  * **minimum_neighbors_generation** (_int_ _|__None_)

  * **minimum_neighbors_value** (_float_ _|__None_)

  * **is_external** (_bool_)

  * **cell_size** (_int_ _|__None_)

  * **metric_values** (_tuple_ _[__float_ _,__...__]__|__None_)




cell_size _: `int` | `None`_ _ = None_#
    

Number of unit cells in the supercell for derivative structure search.

chemical_formula _: `str` | `None`_ _ = None_#
    

Chemical formula of the crystal structure.

_property _chemical_symbols _: list[str] | None_#
    

Return list of chemical symbols in the crystal structure.

fmax _: `float` | `None`_ _ = None_#
    

force convergence criterion in structural optimization.

formation_energy _: `float` | `None`_ _ = None_#
    

Formation energy of the crystal structure in eV (not eV/atom).

_property _formation_energy_per_atom _: float | None_#
    

Return formation energy per atom in eV/atom.

generation _: `int` | `None`_ _ = None_#
    

Generation in genetic algorithm.

get_reduced_composition(_elements_)#
    

Return reduced composition with respect to the given elements.

Parameters:
    

**elements** (`tuple`[`str`, `...`]) – Elements to consider for the reduced composition.

Return type:
    

`list`[`float`] | `None`

Returns:
    

Reduced composition as a list of floats, or None if chemical_formula is not set.

id _: `NewType`(`CrystalStructureID`, `int`)_#
    

Identifier of crystal structure. Expected to be a sequential number.

is_complete()#
    

Return True if all necessary data are set for this crystal structure.

Return type:
    

`bool`

is_external _: `bool`_ _ = False_#
    

True if structure is from external databases like Materials Project.

metric_values _: `tuple`[`float`, `...`] | `None`_ _ = None_#
    

minimum_neighbors_generation _: `int` | `None`_ _ = None_#
    

Used for filtering unpromising crystal structures.

minimum_neighbors_value _: `float` | `None`_ _ = None_#
    

Used for filtering unpromising crystal structures.

_property _mixing_energy _: float | None_#
    

use `formation_energy` instead. Will be removed in MTCSP v2.0.

Type:
    

Deprecated

_property _mixing_energy_per_atom _: float | None_#
    

operator_type _: `str` | `None`_ _ = None_#
    

Type of operator used to generate this crystal structure.

opt_stats _: `dict`[`str`, `float` | `None`] | `None`_ _ = None_#
    

Optimization statistics.

optimized_crystal_structure _: `str` | `None`_ _ = None_#
    

An artifact ID of optimized crystal structure.

optimizer_type _: `Optional`[`Literal`[`'LBFGS'`, `'LBFGSElazy'`, `'FIRE'`]]__ = None_#
    

Type of structural optimizer.

origin _: `str` | `None`_ _ = None_#
    

Descriptor of genetic operation.

parents _: `tuple`[`NewType`(`CrystalStructureID`, `int`), `...`] | `None`_ _ = None_#
    

Identifiers of parent crystal structures.

potential_energy _: `float` | `None`_ _ = None_#
    

Potential energy of the crystal structure with correction terms in eV (not eV/atom).

rejection_reason _: `str` | `None`_ _ = None_#
    

Reason for rejection during local optimization.

steps _: `int` | `None`_ _ = None_#
    

Maximal number of structural optimization steps.

unoptimized_crystal_structure _: `str` | `None`_ _ = None_#
    

An artifact ID of unoptimized crystal structure.

_class _mtcsp.crystal_structure.CrystalStructureID#
    

Identifier of CrystalStructure. Unique within an experiment.

alias of `int`

[ previous mtcsp.convex_hull_search package ](mtcsp.convex_hull_search.html "previous page") [ next mtcsp.derivative_structure package ](mtcsp.derivative_structure.html "next page")

On this page 

  * Module contents
    * `CrystalStructure`
      * `CrystalStructure.cell_size`
      * `CrystalStructure.chemical_formula`
      * `CrystalStructure.chemical_symbols`
      * `CrystalStructure.fmax`
      * `CrystalStructure.formation_energy`
      * `CrystalStructure.formation_energy_per_atom`
      * `CrystalStructure.generation`
      * `CrystalStructure.get_reduced_composition()`
      * `CrystalStructure.id`
      * `CrystalStructure.is_complete()`
      * `CrystalStructure.is_external`
      * `CrystalStructure.metric_values`
      * `CrystalStructure.minimum_neighbors_generation`
      * `CrystalStructure.minimum_neighbors_value`
      * `CrystalStructure.mixing_energy`
      * `CrystalStructure.mixing_energy_per_atom`
      * `CrystalStructure.operator_type`
      * `CrystalStructure.opt_stats`
      * `CrystalStructure.optimized_crystal_structure`
      * `CrystalStructure.optimizer_type`
      * `CrystalStructure.origin`
      * `CrystalStructure.parents`
      * `CrystalStructure.potential_energy`
      * `CrystalStructure.rejection_reason`
      * `CrystalStructure.steps`
      * `CrystalStructure.unoptimized_crystal_structure`
    * `CrystalStructureID`


