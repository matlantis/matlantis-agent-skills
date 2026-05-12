  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.experi...



# mtcsp.experiment package#

## Module contents#

_class _mtcsp.experiment.Experiment(_*_ , _atoms_store_ , _sampler_ , _metric =None_, _relaxer =None_, _experiment_name =None_, _storage =None_, _proxy_host =None_, _proxy_port =None_, _log_file =None_, _rejecter =None_)#
    

Bases: `ExperimentProtocol`

An experiment to search crystal structure.

Parameters:
    

  * **atoms_store** ([_AtomsStore_](mtcsp.atoms.html#mtcsp.atoms.AtomsStore "mtcsp.atoms.AtomsStore"))

  * **sampler** (_CSPSampler_)

  * **metric** (_Metric_ _|__None_)

  * **relaxer** (_Relaxer_ _|__None_)

  * **experiment_name** (_str_ _|__None_)

  * **storage** (_str_ _|__None_)

  * **proxy_host** (_str_ _|__None_)

  * **proxy_port** (_int_ _|__None_)

  * **log_file** (_Path_ _|__None_)

  * **rejecter** (_Rejecter_ _|__None_)




add_crystal_structures_from_experiment(_other_experiment_ , _e_above_hull =0.0_)#
    

Add CrystalStructures from another experiment.

Parameters:
    

  * **other_experiment** (`FrozenExperiment`) – Another instance of Experiment that holds structures to be added.

  * **e_above_hull** (`float`) – Threshold of energy above hull to determine which structures are added. If it’s 0, structures predicted to be thermodynamically stable only are added. (default: `0.0`)



Return type:
    

`list`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]

Returns:
    

Structures that have been added.

_class _mtcsp.experiment.FrozenExperiment(_*_ , _storage_ , _experiment_name_ , _atoms_store =None_, _log_file =None_)#
    

Bases: `ExperimentProtocol`

A frozen experiment.

Parameters:
    

  * **storage** (_str_)

  * **experiment_name** (_str_)

  * **atoms_store** ([_AtomsStore_](mtcsp.atoms.html#mtcsp.atoms.AtomsStore "mtcsp.atoms.AtomsStore") _|__None_)

  * **log_file** (_Path_ _|__None_)




_property _compatibility _: Compatibility | None_#
    

Return the compatibility used in this experiment.

_property _completed_crystal_structures _: list[[CrystalStructure](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]_#
    

Return the list of completed crystal structures.

_property _elements _: Annotated[tuple[str, ...], 'Tuple of chemical symbols']_#
    

Return the elements used in this experiment.

_property _experiment_name _: str_#
    

Return the name of this experiment.

get_atoms_by_id_list(_id_list_)#
    

Retrieves the atoms object of the crystal structure with the specified ID.

Parameters:
    

**id_list** (_list_ _[_[_CrystalStructureID_](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") _]_) – The list of ID of the crystal structure to retrieve.

Returns:
    

The atoms object of the crystal structure with the specified ID.

Return type:
    

[ase.Atoms](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

get_crystal_structure_by_id_list(_id_list_)#
    

Retrieves the crystal structures with the specified ID list.

Parameters:
    

**id_list** (_list_ _[_[_CrystalStructureID_](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") _]_) – The ID list of the crystal structures to retrieve.

Returns:
    

The crystal structures with the specified ID list.

Return type:
    

list[[CrystalStructure](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]

get_crystal_structure_info_by_id_list(_id_list_ , _phase_diagram =None_, _spg_symprec =0.01_, _spg_angle_tolerance =5_)#
    

Retrieves the detailed information of the crystal structure with the specified ID.

Parameters:
    

  * **id_list** (_list_ _[_[_CrystalStructureID_](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructureID "mtcsp.crystal_structure.CrystalStructureID") _]_) – The list of ID of the crystal structure to retrieve.

  * **phase_diagram** ([_MTCSPPhaseDiagram_](mtcsp.analysis.html#mtcsp.analysis.MTCSPPhaseDiagram "mtcsp.analysis.MTCSPPhaseDiagram") _,__optional_) – The phase diagram object to calculate the energy above the hull. Default is None. (default: `None`)

  * **spg_symprec** (_float_ _,__optional_) – Symmetry precision for space group analysis. Default is 0.01. (default: `0.01`)

  * **spg_angle_tolerance** (_float_ _,__optional_) – Angle tolerance for space group analysis. Default is 5. (default: `5`)



Returns:
    

The detailed information of the crystal structure with the specified ID.

Return type:
    

dict[str, str | float]

get_final_stable_crystal_structures()#
    

Return the final stable crystal structures found in this experiment.

Return type:
    

`list`[[`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]

get_final_stable_entries()#
    

Return the final stable entries of `MTCSPEntry` found in this experiment.

Return type:
    

`list`[[`MTCSPEntry`](mtcsp.analysis.html#mtcsp.analysis.MTCSPEntry "mtcsp.analysis._phase_diagram.MTCSPEntry")]

get_state_and_rejection_counts()#
    

Return a dictionary containing counts of trials by their states and rejection reasons.

Return type:
    

`dict`[`str`, `int`]

_property _max_generation _: int_#
    

Return the maximum value of generations.

Note that this might not be the maximum due to the distributed optimization. This might be the maximum minus one.

_property _model_version _: Literal['v0.0.0', 'v2.0.0', 'v3.0.0', 'v4.0.0', 'v5.0.0', 'v6.0.0', 'v7.0.0', 'v8.0.0', 'v9.0.0alpha1+arch19', 'v9.0.0alpha1+arch17', 'vnnp', 'latest'] | None_#
    

Return the model version used in this experiment.

_property _reference_entries _: list[[MTCSPEntry](mtcsp.analysis.html#mtcsp.analysis.MTCSPEntry "mtcsp.analysis._phase_diagram.MTCSPEntry")]_#
    

Return the reference entries used in this experiment.

_class _mtcsp.experiment.ExperimentProtocol(_* args_, _** kwargs_)#
    

Protocol for experiments.

_property _compatibility _: Compatibility | None_#
    

Return the compatibility used in this experiment.

_property _completed_crystal_structures _: list[[CrystalStructure](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")]_#
    

Return the list of completed crystal structures.

_property _elements _: Annotated[tuple[str, ...], 'Tuple of chemical symbols']_#
    

Return the elements of the experiment.

_property _experiment_name _: str_#
    

Return the name of this experiment.

_property _max_generation _: int_#
    

Return the maximum value of generations.

Note that this might not be the maximum due to the distributed optimization. This might be the maximum minus one.

_property _model_version _: Literal['v0.0.0', 'v2.0.0', 'v3.0.0', 'v4.0.0', 'v5.0.0', 'v6.0.0', 'v7.0.0', 'v8.0.0', 'v9.0.0alpha1+arch19', 'v9.0.0alpha1+arch17', 'vnnp', 'latest'] | None_#
    

Return the model version used in this experiment.

_property _reference_entries _: list[[MTCSPEntry](mtcsp.analysis.html#mtcsp.analysis.MTCSPEntry "mtcsp.analysis._phase_diagram.MTCSPEntry")]_#
    

Return the reference entries of the experiment.

[ previous mtcsp.derivative_structure package ](mtcsp.derivative_structure.html "previous page") [ next mtcsp.export package ](mtcsp.export.html "next page")

On this page 

  * Module contents
    * `Experiment`
      * `Experiment.add_crystal_structures_from_experiment()`
    * `FrozenExperiment`
      * `FrozenExperiment.compatibility`
      * `FrozenExperiment.completed_crystal_structures`
      * `FrozenExperiment.elements`
      * `FrozenExperiment.experiment_name`
      * `FrozenExperiment.get_atoms_by_id_list()`
      * `FrozenExperiment.get_crystal_structure_by_id_list()`
      * `FrozenExperiment.get_crystal_structure_info_by_id_list()`
      * `FrozenExperiment.get_final_stable_crystal_structures()`
      * `FrozenExperiment.get_final_stable_entries()`
      * `FrozenExperiment.get_state_and_rejection_counts()`
      * `FrozenExperiment.max_generation`
      * `FrozenExperiment.model_version`
      * `FrozenExperiment.reference_entries`
    * `ExperimentProtocol`
      * `ExperimentProtocol.compatibility`
      * `ExperimentProtocol.completed_crystal_structures`
      * `ExperimentProtocol.elements`
      * `ExperimentProtocol.experiment_name`
      * `ExperimentProtocol.max_generation`
      * `ExperimentProtocol.model_version`
      * `ExperimentProtocol.reference_entries`


  *[*]: Keyword-only parameters separator (PEP 3102)
