  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.analys...



# mtcsp.analysis package#

## Module contents#

### Prototype structure matching#

_class _mtcsp.analysis.MaterialsProjectMatcher(_elements_ , _ltol =0.2_, _stol =0.3_, _angle_tol =5.0_, _mp_api_key =None_)#
    

Bases: `object`

Matches generated structures with structures on the Materials Project.

Parameters:
    

  * **elements** (`Sequence`[`str`]) – List of chemical elements to consider.

  * **ltol** (`float`) – Fractional length tolerance used in [`StructureMatcher`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.structure_matcher.StructureMatcher "\(in pymatgen v2023.7.17\)") (default: `0.2`)

  * **stol** (`float`) – Site tolerance used in [`StructureMatcher`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.structure_matcher.StructureMatcher "\(in pymatgen v2023.7.17\)") (default: `0.3`)

  * **angle_tol** (`float`) – Angle tolerance used in [`StructureMatcher`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.structure_matcher.StructureMatcher "\(in pymatgen v2023.7.17\)") (default: `5.0`)

  * **mp_api_key** (`str` | `None`) – Materials Project API key. If None, will use the default environment variable (default: `None`)




fit_atoms(_atoms_)#
    

Check if the given ASE Atoms matches any Materials Project structure.

Parameters:
    

**atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The ASE Atoms to check.

Return type:
    

[`MPID`](https://materialsproject.github.io/emmet/_reference/emmet.core.mpid.MPID.html#emmet.core.mpid.MPID "\(in emmet\)") | `None`

Returns:
    

The Materials Project ID of the matching structure, or None if no match is found.

fit_structure(_structure_)#
    

Check if the given structure matches any Materials Project structure.

Parameters:
    

**structure** ([`Structure`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.Structure "\(in pymatgen v2023.7.17\)")) – The structure to check.

Return type:
    

[`MPID`](https://materialsproject.github.io/emmet/_reference/emmet.core.mpid.MPID.html#emmet.core.mpid.MPID "\(in emmet\)") | `None`

Returns:
    

The Materials Project ID of the matching structure, or None if no match is found.

### Crystal structure selection#

mtcsp.analysis.select_crystal_structures_with_reduced_composition(_composition_ , _crystal_structure_list_)#
    

Select crystal structures with a specified composition.

Parameters:
    

  * **composition** ([`Composition`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition "\(in pymatgen v2023.7.17\)") | `str`) – Target composition to be matched.

  * **crystal_structure_list** (`Sequence`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]) – List of crystal structures to be filtered



Return type:
    

`list`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]

Returns:
    

List of crystal structures whose compositions are same as the specified composition up to a constant factor.

mtcsp.analysis.select_min_energy_crystal_structures_by_composition(_crystal_structures_)#
    

Select the crystal structure with the lowest formation energy for each composition.

Parameters:
    

**crystal_structures** (`Sequence`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]) – A sequence of crystal structures to be filtered.

Return type:
    

`list`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]

### Phase diagram analysis#

_class _mtcsp.analysis.FormationEnergyEntry(_composition_ , _formation_energy_ , _crystal_structure_id =None_)#
    

Bases: [`PDEntry`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.phase_diagram.PDEntry "\(in pymatgen v2023.7.17\)")

Extension of pymatgen’s PDEntry to represent a formation energy entry.

Be careful self.energy is a formation energy in eV (not eV/atom). For reference simple entries, the formation energy is zero.

Parameters:
    

  * **composition** ([`Composition`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition "\(in pymatgen v2023.7.17\)")) – Composition of the entry.

  * **formation_energy** (`float`) – Formation energy in eV (not eV/atom).

  * **crystal_structure_id** (`Optional`[`NewType`(`CrystalStructureID`, `int`)]) – ID of the crystal structure associated with this entry. (default: `None`)




_property _crystal_structure_id _: [CrystalStructureID](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") | None_#
    

Return the crystal structure ID associated with this entry.

_property _entry_id _: str_#
    

_classmethod _from_crystal_structure(_cs_ , _reference_energies_)#
    

Create FormationEnergyEntry from CrystalStructure and reference energies.

Parameters:
    

  * **cs** ([`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")) – The CrystalStructure to convert.

  * **reference_energies** (`dict`[`str`, `float`]) – A dictionary mapping element symbols to their reference energies (eV/atom).



Return type:
    

`Self`

_classmethod _from_mtcsp_entry(_entry_ , _reference_energies_)#
    

Create FormationEnergyEntry from MTCSPEntry and reference energies.

Parameters:
    

  * **entry** (`MTCSPEntry`) – The MTCSPEntry to convert.

  * **reference_energies** (`dict`[`str`, `float`]) – A dictionary mapping element symbols to their reference energies (eV/atom).



Return type:
    

`Self`

_class _mtcsp.analysis.MixingEnergyEntry(_composition_ , _mixing_energy_ , _crystal_structure_id =None_)#
    

Bases: `FormationEnergyEntry`

Parameters:
    

  * **composition** (_Composition_)

  * **mixing_energy** (_float_)

  * **crystal_structure_id** ([_CrystalStructureID_](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") _|__None_)




_class _mtcsp.analysis.MTCSPEntry(_composition_ , _energy_ , _id =None_)#
    

Bases: [`PDEntry`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.phase_diagram.PDEntry "\(in pymatgen v2023.7.17\)")

Extension of pymatgen’s PDEntry.

Be careful self.energy is not a formation energy but a potential energy in eV.

Parameters:
    

  * **composition** ([`Composition`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.composition.Composition "\(in pymatgen v2023.7.17\)") | `str`) – Composition of the entry.

  * **energy** (`float`) – Potential energy in eV.

  * **id** (`Optional`[`NewType`(`CrystalStructureID`, `int`)]) – ID of the crystal structure associated with this entry. (default: `None`)




as_dict()#
    

Serialize to a dictionary.

Return type:
    

`dict`[`str`, `Any`]

_property _crystal_structure_id _: [CrystalStructureID](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") | None_#
    

Return the crystal structure ID associated with this entry.

_property _entry_id _: str_#
    

_classmethod _from_atoms(_atoms_ , _potential_energy_)#
    

Create MTCSPEntry from ASE Atoms.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The Atoms to convert.

  * **potential_energy** (`float`) – Potential energy in eV.



Return type:
    

`Self`

_classmethod _from_crystal_structure(_cs_)#
    

Create MTCSPEntry from CrystalStructure.

Parameters:
    

**cs** ([`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")) – The CrystalStructure to convert.

Return type:
    

`Self`

_classmethod _from_dict(_dct_)#
    

Deserialize from a dictionary.

Parameters:
    

**dct** (`dict`[`str`, `Any`]) – The dictionary to deserialize.

Return type:
    

`Self`

_class _mtcsp.analysis.MTCSPPhaseDiagram(_entries_ , _*_ , _elements_ , _reference_entries_)#
    

Bases: [`PhaseDiagram`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.phase_diagram.PhaseDiagram "\(in pymatgen v2023.7.17\)")

Extension of pymatgen’s PhaseDiagram class.

Parameters:
    

  * **entries** (`Sequence`[`MTCSPEntry`]) – List of MTCSPEntry to construct the phase diagram.

  * **elements** (`Sequence`[[`Element`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Element "\(in pymatgen v2023.7.17\)") | `str`]) – List of elements in the system.

  * **reference_entries** (`Sequence`[`MTCSPEntry`]) – List of reference entries for each element.




EPS _: `Final`[`float`]__ = 1e-08_#
    

_property _hypervolume _: float_#
    

Return the volume of the energy convex hull below 0 eV region.

_property _reference_energies _: dict[str, float]_#
    

mtcsp.analysis.get_structures_and_e_above_hull(_experiment_ , _max_e_above_hull =None_, _use_reference_entries =False_)#
    

Return list of structures near the convex hull and e_above_hull (eV/atom).

The convex hull is constructed from the experiment.

Parameters:
    

  * **experiment** ([`ExperimentProtocol`](mtcsp.experiment.html#mtcsp.experiment.ExperimentProtocol "mtcsp.experiment._protocol.ExperimentProtocol")) – The experiment to analyze.

  * **max_e_above_hull** (`float` | `None`) – Maximum e_above_hull (eV/atom) to include in the result. (default: `None`)

  * **use_reference_entries** (_bool_)



Return type:
    

`list`[`tuple`[`MTCSPEntry`, `float`]]

Returns:
    

List of tuples of MTCSPEntry and e_above_hull (eV/atom).

### Misc#

mtcsp.analysis.get_crystal_structure_info(_crystal_structure_ , _atoms_ , _phase_diagram =None_, _spg_symprec =0.01_, _spg_angle_tolerance =5.0_)#
    

Extracts and returns detailed information of a given crystal structure.

Parameters:
    

  * **crystal_structure** ([`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")) – The crystal structure to extract information from.

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The atoms object of the crystal structure.

  * **phase_diagram** (`MTCSPPhaseDiagram` | `None`) – The phase diagram object to calculate the energy above the hull. Default is None. (default: `None`)

  * **spg_symprec** (`float`) – Symmetry precision for space group analysis. Default is 0.01. (default: `0.01`)

  * **spg_angle_tolerance** (`float`) – Angle tolerance for space group analysis. Default is 5. (default: `5.0`)



Returns:
    

  * ID

  * Formula

  * Atomic composition fractions of each element

  * Formation energy in eV/atom

  * Lattice parameters (a, b, c, alpha, beta, gamma)

  * Space group number and symbol

  * Crystal system

  * Energy above the hull in eV/atom if phase_diagram is provided




Note that these lattice parameters are for a standardized cell with the given spg_symprec.

Return type:
    

A dictionary containing the extracted crystal structure information, including

mtcsp.analysis.standardize_atoms(_atoms_ , _*_ , _symprec =0.01_, _mag_symprec =None_, _primitive =False_)#
    

Standardize the given Atoms by removing distortions as much as possible.

WARNING: This function may reorder (permute) atoms, so direct index-based comparison might be not valid.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The input atomic structure to standardize

  * **symprec** (`float`) – Symmetry precision threshold for atomic positions (default: 1e-2) (default: `0.01`)

  * **mag_symprec** (`float` | `None`) – Symmetry precision threshold for magnetic moments (default: None, uses 1e-8) (default: `None`)

  * **primitive** (`bool`) – If True, returns primitive cell; if False, returns conventional cell (default: False) (default: `False`)



Return type:
    

[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

Returns:
    

Standardized atoms object with distortions below thresholds removed.

Raises:
    

**ValueError** – If magnetic moments have more than 2 dimensions

mtcsp.analysis.unique_computed_structure_entries(_entries_ , _*_ , _ltol =0.2_, _stol =0.3_, _energy_tol =0.001_)#
    

Unique entries based on geometry and energy.

Parameters:
    

  * **entries** (`list`[[`ComputedStructureEntry`](https://pymatgen.org/pymatgen.entries.html#pymatgen.entries.computed_entries.ComputedStructureEntry "\(in pymatgen v2023.7.17\)")]) – List of computed structure entries to be uniqued.

  * **ltol** (`float`) – Fractional length tolerance used in [`StructureMatcher`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.structure_matcher.StructureMatcher "\(in pymatgen v2023.7.17\)") (default: `0.2`)

  * **stol** (`float`) – Site tolerance used in [`StructureMatcher`](https://pymatgen.org/pymatgen.analysis.html#pymatgen.analysis.structure_matcher.StructureMatcher "\(in pymatgen v2023.7.17\)") (default: `0.3`)

  * **energy_tol** (`float`) – Energy tolerance for uniqueness check, in eV/atom. (default: `0.001`)



Return type:
    

`list`[[`ComputedStructureEntry`](https://pymatgen.org/pymatgen.entries.html#pymatgen.entries.computed_entries.ComputedStructureEntry "\(in pymatgen v2023.7.17\)")]

Returns:
    

List of uniqued computed structure entries.

[ previous Overview ](overview.html "previous page") [ next mtcsp.atoms package ](mtcsp.atoms.html "next page")

On this page 

  * Module contents
    * Prototype structure matching
      * `MaterialsProjectMatcher`
        * `MaterialsProjectMatcher.fit_atoms()`
        * `MaterialsProjectMatcher.fit_structure()`
    * Crystal structure selection
      * `select_crystal_structures_with_reduced_composition()`
      * `select_min_energy_crystal_structures_by_composition()`
    * Phase diagram analysis
      * `FormationEnergyEntry`
        * `FormationEnergyEntry.crystal_structure_id`
        * `FormationEnergyEntry.entry_id`
        * `FormationEnergyEntry.from_crystal_structure()`
        * `FormationEnergyEntry.from_mtcsp_entry()`
      * `MixingEnergyEntry`
      * `MTCSPEntry`
        * `MTCSPEntry.as_dict()`
        * `MTCSPEntry.crystal_structure_id`
        * `MTCSPEntry.entry_id`
        * `MTCSPEntry.from_atoms()`
        * `MTCSPEntry.from_crystal_structure()`
        * `MTCSPEntry.from_dict()`
      * `MTCSPPhaseDiagram`
        * `MTCSPPhaseDiagram.EPS`
        * `MTCSPPhaseDiagram.hypervolume`
        * `MTCSPPhaseDiagram.reference_energies`
      * `get_structures_and_e_above_hull()`
    * Misc
      * `get_crystal_structure_info()`
      * `standardize_atoms()`
      * `unique_computed_structure_entries()`


  *[*]: Keyword-only parameters separator (PEP 3102)
