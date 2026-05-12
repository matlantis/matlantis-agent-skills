# Source code for pfp_api_client.exception
    
    
    import os
    import quopri
    from dataclasses import dataclass
    from datetime import datetime
    from typing import Optional, Sequence, Tuple
    
    
    
    
    [[docs]](../../pfp_api_client.exception.html#pfp_api_client.exception.PFPAPIErrorContext)@dataclass
    class PFPAPIErrorContext:
        code: str
        method: str
        details: str
        exc: Optional[Exception]
        metadata: Sequence[Tuple[str, str]]
        remote: str
        time: datetime
        notebook_id: str
    
        def _qp_value(self, key: str, value: str) -> str:
            if len(key) + len(value) + 2 <= 76:
                # if it fits on one line when adding ": ", fine
                return value
            else:
                # otherwise, assume that it's a long ugly string that we should
                # encode with quoted-printable encoding to avoid copy-paste
                # issues and long lines. make sure that including the "key: "
                # prefix we don't exceed the 76 chars line limit. use quotetabs=True
                # so that we escape tabs and spaces (newlines are preserved).
                dummy_key = "." * (len(key) + 2)
                return quopri.encodestring((dummy_key + value).encode("utf-8"), quotetabs=True).decode(
                    "ascii"
                )[len(dummy_key) :]
    
    
    
    [[docs]](../../pfp_api_client.exception.html#pfp_api_client.exception.PFPAPIErrorContext.format)    def format(self) -> str:
            err = ""
            err += f"time: {self.time.isoformat()}\n"
            err += f"pid: {os.getpid()}\n"
            err += f"code: {self.code}\n"
            err += f"method: {self.method}\n"
            err += f"details: {self._qp_value('details', self.details)}\n"
            err += f"exc: {self._qp_value('exc', repr(self.exc))}\n"
            err += f"notebook_id: {self.notebook_id}\n"
            err += f"metadata: {self._qp_value('metadata', repr(self.metadata))}\n"
            err += f"remote: {self.remote}"
            return err
    
    
    
    
    
    
    [[docs]](../../pfp_api_client.exception.html#pfp_api_client.exception.PFPAPIError)class PFPAPIError(Exception):
        """
        The exception class raised when an error occurs within the PFP API.
        """
    
        def __init__(self, msg: str, context: PFPAPIErrorContext) -> None:
            super().__init__(msg)
            self.msg = msg
            self.context = context
    
        def __str__(self) -> str:
            return (
                self.msg
                + "\n=== When reporting an error, please also share the data below: ===\n"
                + "-" * 76
                + "\n"
                + self.context.format()
                + "\n"
                + "-" * 76
            )
    
    
    
    
    
    
    [[docs]](../../pfp_api_client.exception.html#pfp_api_client.exception.RetriesExceeded)class RetriesExceeded(PFPAPIError):
        """
        The exception class raised when the number of retries specified in the Estimator's \
        `max_retries` is exceeded.
        """
    
        pass
    
    
    
    
    PFCCError = PFPAPIError
    
