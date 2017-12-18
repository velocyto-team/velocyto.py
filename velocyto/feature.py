from typing import *
from collections import defaultdict
import logging
import velocyto as vcy


class Feature:
    """A simple class representing an annotated genomic feature (e.g. exon, intron, masked repeat)"""
    __slots__ = ["start", "end", "kind", "exin_no", "is_validated", "transcript_model"]
    
    def __init__(self, start: int, end: int, kind: int, exin_no: str, transcript_model: Any=None) -> None:
        self.start = start
        self.end = end
        self.transcript_model = transcript_model
        self.kind = kind  # it should be ord("e"), ord("i"), ord("m"), ....
        self.exin_no = int(exin_no)
        self.is_validated = False
    
    def __lt__(self, other: Any) -> bool:
        if self.start == other.start:
            return self.end < other.end
        return self.start < other.start

    def __gt__(self, other: Any) -> bool:
        if self.start == other.start:
            return self.end > other.end
        return self.start > other.start

    def __len__(self) -> int:
        return (self.end - self.start) + 1
        
    def __repr__(self) -> str:
        if self.transcript_model is None:
            return f"Feature not linked to Transcript Model: {self.start}-{self.end} {chr(self.kind)}{self.exin_no}"
        return f"Feature: chr{self.transcript_model.chromstrand}:{self.start}-{self.end} {self.transcript_model.trname}\
    ({self.transcript_model.trid}) {chr(self.kind)}{self.exin_no} {self.transcript_model.genename}({self.transcript_model.geneid})"

    @property
    def is_last_3prime(self) -> bool:
        if self.transcript_model.chromstrand[-1] == "+":
            return self == self.transcript_model.list_features[-1]
        else:
            return self == self.transcript_model.list_features[0]

    def get_downstream_exon(self) -> Any:
        """To use only for introns. Returns the vcy.Feature corresponding to the neighbour exon downstream

        Note
        ----
        In a 15 exons transcript model:
        Downstream to intron10 is exon11 or the interval with index `20` if strand "+".
        Downtream to intron10 is exon10 or the interval with index `10` if strand "-"
        """
        if self.transcript_model.chromstrand[-1] == "+":
            ix = self.exin_no * 2
        else:
            # in the case on strand -
            ix = len(self.transcript_model.list_features) - 2 * self.exin_no + 1
        return self.transcript_model.list_features[ix]

    def get_upstream_exon(self) -> Any:
        """To use only for introns. Returns the vcy.Feature corresponding to the neighbour exon downstream

        Note
        ----
        In a 15 exons transcript model:
        Upstream to intron10 is exon9 or the interval with inxex `18` if strand "+".
        Upstream to intron10 is exon11 or the interval with inxex `8` if strand "-"
        """
        if self.transcript_model.chromstrand[-1] == "+":
            ix = (self.exin_no * 2) - 2
        else:
            # in the case on strand -
            ix = len(self.transcript_model.list_features) - 2 * self.exin_no - 1
        return self.transcript_model.list_features[ix]

    # if self.chromstrand[-1] == "+":
    #             intron_number = self.list_features[-1].exin_no
    #         else:
    #             intron_number = self.list_features[-1].exin_no - 1

    def ends_upstream_of(self, read: vcy.Read) -> bool:
        """The following situation happens
                                                            Read
                                               *|||segment|||-?-||segment|||????????
                ???????|||||Ivl|||||||||*

        """
        return self.end < read.pos  # NOTE: pos is diffetent from start, consider chagning

    def doesnt_start_after(self, segment: Tuple[int, int]) -> bool:
        """One of the following situation happens

                            *||||||segment|||||????????
            *||||Ivl|||||*
                *|||||||||||||Ivl||||||||||????????????
                                    *|||||||||||||Ivl||||||||||????????????
                                              *|||||||||||||Ivl||||||||||????????????

        """
        return self.start < segment[-1]

    def intersects(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        return (segment[-1] - minimum_flanking > self.start) and\
               (segment[0] + minimum_flanking < self.end)  # and ((segment[-1] - segment[0]) > minimum_flanking)

    def contains(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        """One of following situation happens

            *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

                  *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

                      *-----||||||segment|||||-----*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking
        """
        return (segment[0] + minimum_flanking >= self.start) and (segment[-1] - minimum_flanking <= self.end) and ((segment[-1] - segment[0]) > minimum_flanking)

    def start_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        """The following situation happens

          *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking

        """
        return (segment[0] + minimum_flanking < self.start) and (segment[-1] - minimum_flanking > self.start)

    def end_overlaps_with_part_of(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        """The following situation happens

                                      *---|||segment||---*
                *|||||||||||||Ivl||||||||||||||||*

        where `---` idicates the minimum flanking
            
        """
        return (segment[0] + minimum_flanking < self.end) and (segment[-1] - minimum_flanking > self.end)
