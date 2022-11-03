import random
import string
from enum import Enum


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
