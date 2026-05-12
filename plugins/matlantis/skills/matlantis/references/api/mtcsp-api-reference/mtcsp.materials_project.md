  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.materi...



# mtcsp.materials_project package#

## Module contents#

mtcsp.materials_project.download_mp_atoms(_elements_ , _*_ , _max_atoms =None_, _max_e_above_hull =1e-08_, _ignore_simple_elements =True_)#
    

Download Materials Project structures for a given elements.

Parameters:
    

  * **elements** (`Sequence`[`str`]) – Sequence of element symbols to search for.

  * **max_atoms** (`int` | `None`) – Maximum number of atoms in the structure. If None, no limit is applied. (default: `None`)

  * **max_e_above_hull** (`float`) – Maximum energy above hull for the structure. (default: `1e-08`)

  * **ignore_simple_elements** (`bool`) – If True, structures with only one element (e.g., single-element crystals) are ignored (default: `True`)



Return type:
    

`tuple`[`list`[[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")], `list`[`float`], `list`[`str`]]

Returns:
    

List of ASE Atoms objects for the structures. List of formation energies for the structures. List of Materials Project IDs for the structures.

[ previous mtcsp.export package ](mtcsp.export.html "previous page") [ next mtcsp.visualization package ](mtcsp.visualization.html "next page")

On this page 

  * Module contents
    * `download_mp_atoms()`


  *[*]: Keyword-only parameters separator (PEP 3102)
