# Drive#

_class _matlantis_group_drive_client.drive.GroupDriveFileSystem(_* args_, _** kwargs_)[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem)#
    

Bases: `AbstractFileSystem`

DEFAULT_CHUNK_SIZE _ = 5242880_#
    

DELIMITER _ = '/'_#
    

FILTERED_KEYS _ = ('created', 'drive_name', 'last_modified', 'name', 'type', 'size')_#
    

cat_file(_path : str_, _start : Optional[int] = None_, _end : Optional[int] = None_, _** kwargs: Kwargs_) → bytes[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.cat_file)#
    

Get the content of a file

Parameters:
    

  * **path** (_URL of file on this filesystems_) – 

  * **start** (_int_) – Bytes limits of the read. If negative, backwards from end, like usual python slices. Either can be None for start or end of file, respectively

  * **end** (_int_) – Bytes limits of the read. If negative, backwards from end, like usual python slices. Either can be None for start or end of file, respectively

  * **kwargs** (passed to `open()`.) – 




download(_drive_path : str_, _local_path : str_) → None[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.download)#
    

Alias of AbstractFileSystem.get.

drive_list() → List[Dict[str, str]][[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.drive_list)#
    

ls(_path : str_, _detail : bool = True_, _** kwargs: Kwargs_) → List[Dict[str, str]][[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.ls)#
    

List objects at path.

This should include subdirectories and files at that location. The difference between a file and a directory must be clear when details are requested.

The specific keys, or perhaps a FileInfo class, or similar, is TBD, but must be consistent across implementations. Must include:

  * full path to the entry (without protocol)

  * size of the entry, in bytes. If the value cannot be determined, will be `None`.

  * type of entry, “file”, “directory” or other




Additional information may be present, aproriate to the file-system, e.g., generation, checksum, etc.

May use refresh=True|False to allow use of self._ls_from_cache to check for a saved listing and avoid calling the backend. This would be common where listing may be expensive.

Parameters:
    

  * **path** (_str_) – 

  * **detail** (_bool_) – if True, gives a list of dictionaries, where each is the same as the result of `info(path)`. If False, gives a list of paths (str).

  * **kwargs** (_may have additional backend-specific options_ _,__such as version_) – information



Returns:
    

  * _List of strings if detail is False, or list of directory information_

  * _dicts if detail is True._




mkdir(_path : str_, _create_parents : bool = True_, _** kwargs: Kwargs_) → None[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.mkdir)#
    

Create directory entry at path

For systems that don’t have true directories, may create an for this instance only and not touch the real filesystem

Parameters:
    

  * **path** (_str_) – location

  * **create_parents** (_bool_) – if True, this is equivalent to `makedirs`

  * **kwargs** – may be permissions, etc.




pipe_file(_path : str_, _value : bytes_, _** _: Kwargs_) → None[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.pipe_file)#
    

Set the bytes of given file

upload(_local_path : str_, _drive_path : str_, _chunk_size : int = 5242880_) → None[[source]](../_modules/matlantis_group_drive_client/drive.html#GroupDriveFileSystem.upload)#
    

Alias of AbstractFileSystem.put.
