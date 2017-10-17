from typing import *
import inspect


class Read:
    """ Container for reads from sam alignment file """
    __slots__ = ["bc", "umi", "chrom", "strand", "pos", "cigar", "segments", "clip5", "clip3", "spliced"]

    def __init__(self, bc: str, umi: str, chrom: str, strand: str, pos: int, cigar: str, segments: List, clip5: Any, clip3: Any, spliced: bool) -> None:
        self.bc, self.umi, self.chrom, self.strand, self.pos, self.cigar, self.segments, self.clip5, self.clip3, self.spliced = \
            bc, umi, chrom, strand, pos, cigar, segments, clip5, clip3, spliced

    def is_spliced(self) -> bool:
        return self.spliced  # len(self.segments) > 1
        
    def start(self) -> int:
        return self.segments[0][0]

    def end(self) -> int:
        return self.segments[-1][1]

    def __str__(self) -> str:
        tmp = ""
        for i in self.__slots__:
            attribute = getattr(self, i)
            tmp += f"{i}: {type(attribute).__name__}={attribute}, "
        return tmp
