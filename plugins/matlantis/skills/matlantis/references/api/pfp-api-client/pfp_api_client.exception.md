# Exceptions#

Exception classes are provided for PFP-related errors.

pfp_api_client.exception.PFCCError#
    

alias of `PFPAPIError`

_exception _pfp_api_client.exception.PFPAPIError(_msg : str_, _context : PFPAPIErrorContext_)[[source]](_modules/pfp_api_client/exception.html#PFPAPIError)#
    

Bases: `Exception`

The exception class raised when an error occurs within the PFP API.

_class _pfp_api_client.exception.PFPAPIErrorContext(_code : str_, _method : str_, _details : str_, _exc : Optional[Exception]_, _metadata : Sequence[Tuple[str, str]]_, _remote : str_, _time : datetime.datetime_, _notebook_id : str_)[[source]](_modules/pfp_api_client/exception.html#PFPAPIErrorContext)#
    

Bases: `object`

code _: str_#
    

details _: str_#
    

exc _: Optional[Exception]_#
    

format() → str[[source]](_modules/pfp_api_client/exception.html#PFPAPIErrorContext.format)#
    

metadata _: Sequence[Tuple[str, str]]_#
    

method _: str_#
    

notebook_id _: str_#
    

remote _: str_#
    

time _: datetime_#
    

_exception _pfp_api_client.exception.RetriesExceeded(_msg : str_, _context : PFPAPIErrorContext_)[[source]](_modules/pfp_api_client/exception.html#RetriesExceeded)#
    

Bases: `PFPAPIError`

The exception class raised when the number of retries specified in the Estimator’s max_retries is exceeded.

[ __ previous Utilities ](pfp_api_client.pfp.utils.html "previous page")
