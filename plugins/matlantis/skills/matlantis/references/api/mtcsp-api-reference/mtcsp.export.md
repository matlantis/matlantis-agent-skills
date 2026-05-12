  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.export package



# mtcsp.export package#

## Module contents#

mtcsp.export.export_as_computed_structure_entry(_cs_ , _atoms_store_)#
    

Export crystal structure as pymatgen’s [`ComputedStructureEntry`](https://pymatgen.org/pymatgen.entries.html#pymatgen.entries.computed_entries.ComputedStructureEntry "\(in pymatgen v2023.7.17\)").

Parameters:
    

  * **cs** ([`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure")) – Crystal structure to export.

  * **atoms_store** ([`AtomsStore`](mtcsp.atoms.html#mtcsp.atoms.AtomsStore "mtcsp.atoms.AtomsStore")) – Atoms store to load the optimized crystal structure.



Returns:
    

The entry’s `data` dict contains `"crystal_structure_id"` (`cs.id`) for mapping back to the original [`CrystalStructure`](mtcsp.crystal_structure.html#mtcsp.crystal_structure.CrystalStructure "mtcsp.crystal_structure.CrystalStructure"). Returns None when the atoms are missing in the store.

Return type:
    

[`ComputedStructureEntry`](https://pymatgen.org/pymatgen.entries.html#pymatgen.entries.computed_entries.ComputedStructureEntry "\(in pymatgen v2023.7.17\)") | None

mtcsp.export.export_experiment(_experiment_name_ , _storage_ , _atoms_store_)#
    

Export the experiment to separate sqlite file and export the atoms in the experiment to a zip file.

The following two files are created:

  * {experiment_name}.sqlite: The sqlite file of the experiment.

  * {experiment_name}.zip: The zip file of the atoms in the experiment.




Parameters:
    

  * **experiment_name** (_str_) – Name of the experiment to export.

  * **storage** (`StorageLikeStr`) – Storage of the experiment to export.

  * **atoms_store** (`AtomsStore`)



Return type:
    

`None`

[ previous mtcsp.experiment package ](mtcsp.experiment.html "previous page") [ next mtcsp.materials_project package ](mtcsp.materials_project.html "next page")

On this page 

  * Module contents
    * `export_as_computed_structure_entry()`
    * `export_experiment()`


