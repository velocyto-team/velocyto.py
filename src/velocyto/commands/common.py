import random
import string
from datetime import datetime
from enum import Enum
from sys import stderr

from loguru import logger

from velocyto import logic


class logicType(str, Enum):
    Permissive10X = "Permissive10X"
    Intermediate10X = "Intermediate10X"
    ValidatedIntrons10X = "ValidatedIntrons10X"
    Stricter10X = "Stricter10X"
    ObservedSpanning10X = "ObservedSpanning10X"
    Discordant10X = "Discordant10X"
    SmartSeq2 = "SmartSeq2"


class UMIExtension(str, Enum):
    no = "no"
    char = "chr"
    gene = "Gene"
    Nbp = "[N]bp"


class loomdtype(str, Enum):
    uint16 = "uint16"
    uint32 = "uint32"
    uint64 = "uint64"


def id_generator(size: int = 6, chars: str = string.ascii_uppercase + string.digits) -> str:
    return "".join(random.choice(chars) for _ in range(size))


def choose_logic(choice: logicType) -> logic.Logic:
    if choice == "Permissive10X":
        return logic.Permissive10X
    elif choice == "Intermediate10X":
        return logic.Intermediate10X
    elif choice == "ValidatedIntrons10X":
        return logic.ValidatedIntrons10X
    elif choice == "Stricter10X":
        return logic.Stricter10X
    elif choice == "ObservedSpanning10X":
        return logic.ObservedSpanning10X
    elif choice == "Discordant10X":
        return logic.Discordant10X
    elif choice == "SmartSeq2":
        return logic.SmartSeq2
    else:
        logger.error(f"{choice.value} is not a valid logic type")
        exit()


def choose_dtype(choice: loomdtype) -> str:
    if choice == "uint16":
        return "uint16"
    elif choice == "uint32":
        return "uint32"
    else:
        return "uint64"


def init_logger(verbose: int, msg_format: str = None) -> None:
    if msg_format is None:
        msg_format = "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>·-·<level>{message}</level>"

    logger.add(f"velocyto_{datetime.now().strftime('%d-%m-%Y--%H-%M-%S')}.log", level="DEBUG")

    if verbose == 3:
        logger.add(stderr, format=msg_format, level="DEBUG")
    elif verbose == 2:
        logger.add(stderr, format=msg_format, level="INFO")
    elif verbose == 1:
        logger.add(stderr, format=msg_format, level="WARNING")
    else:
        logger.add(stderr, format=msg_format, level="ERROR")
