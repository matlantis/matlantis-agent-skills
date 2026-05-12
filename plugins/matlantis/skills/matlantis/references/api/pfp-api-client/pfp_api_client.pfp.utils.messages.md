# Messages#

_class _pfp_api_client.pfp.utils.messages.MessageEnum(_value_ , _names =None_, _*_ , _module =None_, _qualname =None_, _type =None_, _start =1_, _boundary =None_)[[source]](_modules/pfp_api_client/pfp/utils/messages.html#MessageEnum)#
    

Bases: `IntEnum`

Enum class which defines the messages used in Estimator.

Variables
    

  * **ExperimentalElementWarning** – Experimental element was detected.

  * **UnexpectedElementWarning** – Unexpected element was detected.




ExperimentalElementWarning _ = 1001_#
    

UnexpectedElementWarning _ = 1002_#
    

_class _pfp_api_client.pfp.utils.messages.MessageType(_value_ , _names =None_, _*_ , _module =None_, _qualname =None_, _type =None_, _start =1_, _boundary =None_)[[source]](_modules/pfp_api_client/pfp/utils/messages.html#MessageType)#
    

Bases: `IntEnum`

Enum class which defines the message type.

CustomMessage _ = 10_#
    

CustomWarning _ = 20_#
    

pfp_api_client.pfp.utils.messages.message_info(_message : MessageEnum_) → Tuple[MessageType, str][[source]](_modules/pfp_api_client/pfp/utils/messages.html#message_info)#
    

Return the message type and corresponding string.

Parameters
    

**message** (_MessageEnum_) – Input message.

Returns
    

**results**

Return type
    

Tuple[MessageType, str]
