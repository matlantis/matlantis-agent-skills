  * [ ](index.html)
  * [MTCSP Reference](mtcsp.html)
  * mtcsp.atoms package



# mtcsp.atoms package#

## Module contents#

_class _mtcsp.atoms.AtomsStore(_* args_, _** kwargs_)#
    

A protocol defining the interface for an atoms storage backend.

load(_atoms_id_)#
    

Load and return the atoms corresponding to the given artifact ID.

Parameters:
    

**atoms_id** (`str`) – The artifact ID of the atoms to be loaded.

Return type:
    

[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

Returns:
    

The loaded ASE Atoms object.

save(_atoms_ , _artifact_id =None_)#
    

Save the given atoms and return its artifact ID.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The ASE Atoms object to be saved.

  * **artifact_id** (`str` | `None`) – The artifact ID to use for saving the atoms. If None, a new UUID-based ID will be generated. (default: `None`)



Return type:
    

`str`

Returns:
    

The artifact ID under which the atoms were saved.

_class _mtcsp.atoms.AtomsStoreWithBackend(_backend_)#
    

Bases: `AtomsStore`

Parameters:
    

**backend** (_ArtifactStore_)

load(_atoms_id_)#
    

Load and return the atoms corresponding to the given artifact ID.

Parameters:
    

**atoms_id** (`str`) – The artifact ID of the atoms to be loaded.

Return type:
    

[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

Returns:
    

The loaded ASE Atoms object.

save(_atoms_ , _artifact_id =None_)#
    

Save the given atoms and return its artifact ID.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The ASE Atoms object to be saved.

  * **artifact_id** (`str` | `None`) – The artifact ID to use for saving the atoms. If None, a new UUID-based ID will be generated. (default: `None`)



Return type:
    

`str`

Returns:
    

The artifact ID under which the atoms were saved.

_class _mtcsp.atoms.FileSystemAtomsStore(_base_path_)#
    

Bases: `AtomsStoreWithBackend`

Atoms store that saves atoms as JSON files in a filesystem directory.

If the number of atoms files is too large, consider using `ZipFilesAtomsStore` instead.

Parameters:
    

**base_path** (`str` | `Path`) – The base directory where atoms JSON files will be stored.

_class _mtcsp.atoms.ZipFileAtomsStore(_zip_path_ , _safemode =False_)#
    

Bases: `AtomsStore`

Atoms store that saves atoms as JSON files in a single Zip file.

Parameters:
    

  * **zip_path** (`str` | `Path`) – The path to the Zip file where atoms JSON files will be stored.

  * **safemode** (`bool`) – If True, avoid directly appending contents to the Zip file to prevent corruption. (default: `False`)




load(_atoms_id_)#
    

Load and return the atoms corresponding to the given artifact ID.

Parameters:
    

**atoms_id** (`str`) – The artifact ID of the atoms to be loaded.

Return type:
    

[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

Returns:
    

The loaded ASE Atoms object.

save(_atoms_ , _artifact_id =None_)#
    

Save the given atoms and return its artifact ID.

Parameters:
    

  * **atoms** ([`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")) – The ASE Atoms object to be saved.

  * **artifact_id** (`str` | `None`) – The artifact ID to use for saving the atoms. If None, a new UUID-based ID will be generated. (default: `None`)



Return type:
    

`str`

Returns:
    

The artifact ID under which the atoms were saved.

_class _mtcsp.atoms.ZipFilesAtomsStore(_directory_)#
    

Bases: `AtomsStore`

Writes atoms json files into a Zip file per process.

Please use this class when you have a large number of atoms files to store.

Parameters:
    

**directory** (`str` | `Path`) – The directory where zip files are stored.

load(_atoms_id_)#
    

Loads atoms from zip files in the directory.

Tries to find {atoms_id}.json in all zip files sequentially.

Return type:
    

[`Atoms`](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)")

Parameters:
    

**atoms_id** (_str_)

save(_atoms_ , _artifact_id =None_)#
    

Append an atoms json file to the zip file dedicated to the process.

Return type:
    

`str`

Parameters:
    

  * **atoms** ([_Atoms_](https://ase-lib.org/ase/atoms.html#ase.Atoms "\(in ASE\)"))

  * **artifact_id** (_str_ _|__None_)




[ previous mtcsp.analysis package ](mtcsp.analysis.html "previous page") [ next mtcsp.convex_hull_search package ](mtcsp.convex_hull_search.html "next page")

On this page 

  * Module contents
    * `AtomsStore`
      * `AtomsStore.load()`
      * `AtomsStore.save()`
    * `AtomsStoreWithBackend`
      * `AtomsStoreWithBackend.load()`
      * `AtomsStoreWithBackend.save()`
    * `FileSystemAtomsStore`
    * `ZipFileAtomsStore`
      * `ZipFileAtomsStore.load()`
      * `ZipFileAtomsStore.save()`
    * `ZipFilesAtomsStore`
      * `ZipFilesAtomsStore.load()`
      * `ZipFilesAtomsStore.save()`


