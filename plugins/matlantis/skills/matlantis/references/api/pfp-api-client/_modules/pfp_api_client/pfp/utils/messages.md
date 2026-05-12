# Source code for pfp_api_client.pfp.utils.messages
    
    
    from enum import IntEnum
    from typing import Tuple
    
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageEnum)class MessageEnum(IntEnum):
        """
        Enum class which defines the messages used in `Estimator`.
    
        :ivar ExperimentalElementWarning: Experimental element was detected.
        :ivar UnexpectedElementWarning: Unexpected element was detected.
        """
    
        ExperimentalElementWarning = 1001
        UnexpectedElementWarning = 1002
    
        def __str__(self) -> str:
            message_type, message_str = message_info(self)
            if message_type == MessageType.CustomMessage:
                return "Runtime message (ID {:d}): {:s}".format(int(self), message_str)
            elif message_type == MessageType.CustomWarning:
                return "Runtime warning (ID {:d}): {:s}".format(int(self), message_str)
            else:
                raise ValueError("Unknown MessageType")
    
    
    
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.MessageType)class MessageType(IntEnum):
        """
        Enum class which defines the message type.
        """
    
        CustomWarning = 20
        CustomMessage = 10
    
    
    
    
    _CW = MessageType.CustomWarning
    _CM = MessageType.CustomMessage
    _MESSAGE_DICT = {
        MessageEnum.ExperimentalElementWarning: (
            _CW,
            "Experimental element was detected.\
             You can suppress this message with Estimator.set_message_status().",
        ),
        MessageEnum.UnexpectedElementWarning: (
            _CW,
            "Unexpected element was detected. There might be strange behavior.",
        ),
    }
    
    
    
    
    [[docs]](../../../../pfp_api_client.pfp.utils.messages.html#pfp_api_client.pfp.utils.messages.message_info)def message_info(message: MessageEnum) -> Tuple[MessageType, str]:
        """
        Return the message type and corresponding string.
    
        Parameters
        --------
        message: MessageEnum
            Input message.
    
        Returns
        --------
        results: Tuple[MessageType, str]
        """
        return _MESSAGE_DICT[message]
    
    
    
