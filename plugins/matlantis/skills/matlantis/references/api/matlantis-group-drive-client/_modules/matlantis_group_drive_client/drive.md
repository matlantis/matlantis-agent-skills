# Source code for matlantis_group_drive_client.drive
    
    
    import base64
    import dataclasses
    import io
    import json
    import os
    import time
    from typing import Any, Dict, Iterator, List, NewType, Optional, SupportsBytes, Union, cast
    from urllib.parse import urljoin
    
    import requests
    from dacite.core import from_dict
    from filelock import BaseFileLock, FileLock
    from fsspec.spec import AbstractBufferedFile, AbstractFileSystem
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    
    from matlantis_group_drive_client.exception import DirectoryError
    
    Kwargs = NewType("Kwargs", Dict[str, Any])
    
    
    @dataclasses.dataclass
    class _Token:
        access_token: str
        # client does not need refresh_token currently
        # refresh_token: str
        expires_at: int
    
    
    @dataclasses.dataclass
    class _Config:
        url: str
        token: Optional[_Token] = None
        k8s_jwt: Optional[str] = None
    
    
    class _ConfigSession(requests.Session):
        """Wrapper class of requests.Session
        It reads config file and updates access_token regularly when it requests (via self.get(...), self.post(...), etc.).
        """
    
        CONFIG_PATH_ENV_NAME = "MATLANTIS_CONTENTS_PRIVATE_API_CONFIG_PATH"
        LOAD_CONFIG_INTERVAL = 5 * 60  # 5 minutes
        GROUP_DRIVE_CLIENT_AUTH_V2_ENV_NAME = "GROUP_DRIVE_CLIENT_AUTH_V2"
        MATLANTIS_CONTENTS_PRIVATE_API_URL_ENV_NAME = "MATLANTIS_CONTENTS_PRIVATE_API_URL"
        K8S_JWT_PATH = "/var/run/secrets/tokens/notebook-user-token"
        DELIMITER = "/"
    
        _lock: BaseFileLock
        _config: Optional[_Config]
        _last_updated: float
    
        def __init__(self, lock_timeout: float = 1) -> None:
            super().__init__()
            self._config = None
            if self._is_auth_v2_enabled():
                self._load_k8s_jwt_if_necessary()
            else:
                lock_file_path = self.config_path + ".lock"
                self._lock = FileLock(lock_file_path, timeout=lock_timeout)
                self._load_config_if_necessary()
    
            # Setup retry strategy
            retry_strategy = Retry(
                total=3,
                backoff_factor=1,
                status_forcelist=[429, 500, 502, 503, 504],
                allowed_methods=["HEAD", "GET", "POST", "PUT", "DELETE", "OPTIONS", "TRACE"],
            )
            adapter = HTTPAdapter(max_retries=retry_strategy)
            self.mount("http://", adapter)
            self.mount("https://", adapter)
    
        def _load_config_if_necessary(self) -> None:
            if self._config is not None and time.time() - self._last_updated < self.LOAD_CONFIG_INTERVAL:
                return
            with self._lock:
                with open(self.config_path, "r") as f:
                    self._config = from_dict(data_class=_Config, data=json.load(f))
                    self._last_updated = time.time()
    
        def _load_k8s_jwt_if_necessary(self) -> None:
            if self._config is not None and time.time() - self._last_updated < self.LOAD_CONFIG_INTERVAL:
                return
            with open(self.K8S_JWT_PATH) as f:
                self._config = _Config(
                    url=os.getenv(self.MATLANTIS_CONTENTS_PRIVATE_API_URL_ENV_NAME), k8s_jwt=f.read().strip()
                )
                self._last_updated = time.time()
    
        def _is_auth_v2_enabled(self) -> bool:
            return os.getenv(self.GROUP_DRIVE_CLIENT_AUTH_V2_ENV_NAME) == "1"
    
        @property
        def config_path(self) -> str:
            _path = os.getenv(self.CONFIG_PATH_ENV_NAME)
            if _path is None:
                raise ValueError("Config path is invalid.")
            return _path
    
        @property
        def config(self) -> _Config:
            if self._is_auth_v2_enabled():
                self._load_k8s_jwt_if_necessary()
            else:
                self._load_config_if_necessary()
            assert self._config is not None
            return self._config
    
        @property
        def url(self) -> str:
            return self.config.url.rstrip(self.DELIMITER) + self.DELIMITER
    
        def request(  # type: ignore[override]
            self,
            method: str,
            url: Union[str, bytes],
            *args: Any,
            **kwargs: Any,
        ) -> requests.Response:
            if self._is_auth_v2_enabled():
                self._load_k8s_jwt_if_necessary()
                self.headers["Authorization"] = f"Bearer {self.config.k8s_jwt}"
            else:
                self._load_config_if_necessary()
                self.headers["Authorization"] = f"Bearer {self.config.token.access_token}"
            full_url: Union[str, bytes]
            if isinstance(url, str):
                full_url = urljoin(self.url, url)
            else:
                full_url = urljoin(self.url.encode(), url)
            return super().request(method, full_url, *args, **kwargs)
    
    
    class _GroupDriveFileWriter(AbstractBufferedFile):  # type:ignore[misc]
        FILES_URL = "files"
        MULTIPART_URL = "files/multipart"
    
        _chunk: int  # chunk should be from 1 to 10000
        _config_session: _ConfigSession
        _upload_id: Optional[str]
    
        def __init__(
            self, fs: AbstractFileSystem, path: str, session: _ConfigSession, mode: str = "wb", **kwargs: Kwargs
        ) -> None:
            super().__init__(fs=fs, path=path, mode=mode, cache_type="none", **kwargs)
            self._chunk = 1
            self._config_session = session
            self._upload_id = None
            if mode != "wb":
                raise ValueError(f"Invalid mode: {mode}")
    
        def _initiate_upload(self) -> None:
            response = self._config_session.post(self.MULTIPART_URL, params=dict(content_path=self.path))
            response.raise_for_status()
    
            upload_id = response.json().get("upload_id")
            if not upload_id:
                raise ValueError("upload_id is invalid")
            self._upload_id = upload_id
    
        def _upload_chunk(self, final: bool = False) -> bool:
            if final:
                # Need to set -1 to chunk when this is a final chunk
                chunk = -1
            else:
                chunk = self._chunk
                self._chunk += 1
    
            self.buffer.seek(0)
            data = self.buffer.read()
            response = self._config_session.put(
                self.FILES_URL, params=dict(content_path=self.path, chunk=chunk, upload_id=self._upload_id), data=data
            )
            response.raise_for_status()
            return True
    
    
    class _GroupDriveStreamFile(AbstractBufferedFile):  # type:ignore[misc]
        DEFAULT_BLOCK_SIZE = 5 * 1024 * 1024  # 5MB
        FILES_URL = "files"
    
        fs: "GroupDriveFileSystem"
        _config_session: _ConfigSession
        _stream: requests.Response
        _iter_content: Iterator[SupportsBytes]
    
        def __init__(
            self, fs: "GroupDriveFileSystem", path: str, session: _ConfigSession, mode: str = "rb", **kwargs: Kwargs
        ) -> None:
            super().__init__(fs=fs, path=path, mode=mode, cache_type="none", **kwargs)
            self._config_session = session
            if mode != "rb":
                raise ValueError(f"Invalid mode: {mode}")
            self._stream = self._config_session.get(self.FILES_URL, params=dict(content_path=path), stream=True)
            self._stream.raise_for_status()
            self._iter_content = self._stream.iter_content(chunk_size=self.DEFAULT_BLOCK_SIZE)
    
        def seek(self, offset: int, *args: Any, **kwargs: Kwargs) -> None:
            # Allow to seek if the offset equals to tell()
            if offset == self.tell():
                return
            raise ValueError("Cannot seek streaming file on Group Drive")
    
        def read(self, num: int = -1) -> bytes:
            if num != -1:
                # This does not support to specify reading position
                # because Group Drive does not support file seek
                raise ValueError("Group Drive does not support file seek so you cannot specify reading position")
            return self.fs.cat_file(self.path)
    
        def __next__(self) -> bytes:
            return bytes(self._iter_content.__next__())
    
        def close(self) -> None:
            self._stream.close()
            super().close()
    
        def __eq__(self, other: Any) -> bool:
            # FIXME: More elaborate implementation
            if "mode" not in other:
                return False
            return cast(bool, super().__eq__(other))
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem)class GroupDriveFileSystem(AbstractFileSystem):  # type:ignore[misc]
        DEFAULT_CHUNK_SIZE = 5 * 1024 * 1024  # 5MB
        FILTERED_KEYS = ("created", "drive_name", "last_modified", "name", "type", "size")
        DELIMITER = "/"
    
        _config_session: _ConfigSession
    
        def __init__(self, lock_timeout: float = 1, **kwargs: Kwargs) -> None:
            self._config_session = _ConfigSession(lock_timeout=lock_timeout)
            super().__init__(**kwargs)
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.ls)    def ls(self, path: str, detail: bool = True, **kwargs: Kwargs) -> List[Dict[str, str]]:
            return [self._filter_keys(row) for row in self._ls(path=path, **kwargs)]
    
    
    
        def _ls(self, path: str, **kwargs: Kwargs) -> List[Dict[str, str]]:
            content = "1"
            r = self._config_session.get("", params=dict(content_path=path, content=content))
            r.raise_for_status()
    
            out = r.json()
    
            if out["type"] == "directory":
                out = out["content"]
            else:
                out = [out]
    
            if out is None:
                return []
    
            for o in out:
                o["name"] = o.pop("path")
                o.pop("content")
                if o["type"] == "notebook":
                    o["type"] = "file"
            return cast(List[Dict[str, str]], out)
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.cat_file)    def cat_file(self, path: str, start: Optional[int] = None, end: Optional[int] = None, **kwargs: Kwargs) -> bytes:
            r = self._config_session.get("", params=dict(content_path=path, content="1"))
            r.raise_for_status()
            out = r.json()
    
            if out["type"] == "directory":
                raise DirectoryError(path)
    
            if out["format"] == "text":
                b = str(out["content"]).encode()
            elif out["format"] == "json":
                b = json.dumps(out["content"]).encode()
            else:
                b = base64.b64decode(out["content"])
            return b[start:end]
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.pipe_file)    def pipe_file(self, path: str, value: bytes, **_: Kwargs) -> None:
            path = self._strip_protocol(path)
            json = {
                "name": path.rsplit(self.DELIMITER, 1)[-1],
                "path": path,
                "size": len(value),
                "content": base64.b64encode(value).decode(),
                "format": "base64",
                "type": "file",
            }
            r = self._config_session.put("", params=dict(content_path=path), json=json)
            r.raise_for_status()
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.mkdir)    def mkdir(self, path: str, create_parents: bool = True, **kwargs: Kwargs) -> None:
            if create_parents and self.DELIMITER in path:
                self.mkdir(path.rsplit(self.DELIMITER, 1)[0], True)
            json = {
                "name": path.rsplit(self.DELIMITER, 1)[-1],
                "path": path,
                "size": 0,
                "content": "",
                "format": "",
                "type": "directory",
            }
            r = self._config_session.put("", params=dict(content_path=path), json=json)
            r.raise_for_status()
    
    
    
        def _rm(self, path: str) -> None:
            response = self._config_session.delete("", params=dict(content_path=path))
            response.raise_for_status()
    
        def _open(self, path: str, mode: str = "rb", **kwargs: Kwargs) -> Union[io.BytesIO, _GroupDriveFileWriter]:
            if mode == "rb":
                return _GroupDriveStreamFile(self, path=path, session=self._config_session, mode=mode)
            elif mode == "wb":
                return _GroupDriveFileWriter(self, path=path, session=self._config_session, mode=mode)
            else:
                raise ValueError(f"Invalid mode: {mode} (must be 'rb' or 'wb')")
    
        def _filter_keys(self, row: Dict[str, str]) -> Dict[str, str]:
            return {key: row[key] for key in row.keys() & self.FILTERED_KEYS}
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.drive_list)    def drive_list(self) -> List[Dict[str, str]]:
            return self.ls("")
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.upload)    def upload(self, local_path: str, drive_path: str, chunk_size: int = DEFAULT_CHUNK_SIZE) -> None:
            with open(local_path, "rb") as f_local:
                with self.open(drive_path, "wb") as f_drive:
                    while True:
                        chunk_data = f_local.read(chunk_size)
                        f_drive.write(chunk_data)
                        if len(chunk_data) == 0:
                            return
    
    
    
    
    
    [[docs]](../../modules/matlantis_group_drive_client.drive.html#matlantis_group_drive_client.drive.GroupDriveFileSystem.download)    def download(self, drive_path: str, local_path: str) -> None:
            with self.open(drive_path, "rb") as f_drive:
                with open(local_path, "wb") as f_local:
                    for data in f_drive:
                        f_local.write(data)
    
    
    
