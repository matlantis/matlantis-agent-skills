# Source code for matlantis_group_drive_client.exception
    
    
    from typing import Optional
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.FileError)class FileError(Exception):
        message = "File not found"
    
        def __init__(self, path: Optional[str] = None) -> None:
            if path is not None:
                self.message = f"{self.message}: {path}"
    
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.FileNotFoundError)class FileNotFoundError(FileError):
        message = "File not found"
    
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.exception.html#matlantis_group_drive_client.exception.DirectoryError)class DirectoryError(FileError):
        message = "This is a directory"
    
    
    
