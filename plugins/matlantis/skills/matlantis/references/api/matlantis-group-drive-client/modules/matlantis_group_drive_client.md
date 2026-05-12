# Matlantis Group Drive Client Reference#

This library is implemented based on a file system access interface called fsspec. See also [the documentation of fsspec](https://filesystem-spec.readthedocs.io/en/latest/).

Note

Group Drive is inherently different from normal file systems such as NFS in order to provide higher availability, scalability and other advantages. On the other hands, some functions such as file seek cannot be used.

## Drive#

This module contains the core file system object for Matlantis Group Drive.

[`GroupDriveFileSystem`](matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem "matlantis_group_drive_client.drive.GroupDriveFileSystem")(*args, **kwargs) |   
---|---  
  
## Exception#

This module contains the exception classes

[`FileError`](matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.FileError "matlantis_group_drive_client.exception.FileError")([path]) |   
---|---  
[`FileNotFoundError`](matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.FileNotFoundError "matlantis_group_drive_client.exception.FileNotFoundError")([path]) |   
[`DirectoryError`](matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.DirectoryError "matlantis_group_drive_client.exception.DirectoryError")([path]) |   
  
[ __ previous Welcome to Matlantis Group Drive Client documentation! ](../index.html "previous page")

__On this page

  * Drive 
  * Exception 


