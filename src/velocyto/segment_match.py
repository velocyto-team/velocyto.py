from .constants import SPLIC_INACUR
from .feature import Feature


class SegmentMatch:
    __slots__ = ["segment", "feature", "is_spliced"]

    def __init__(self, segment: tuple[int, int], feature: Feature, is_spliced: bool = False) -> None:
        self.segment = segment
        self.feature = feature
        self.is_spliced = is_spliced  # this is really BAM_CREF_SKIP

    @property
    def maps_to_intron(self) -> bool:
        return self.feature.kind == 105  # ord("i")

    @property
    def maps_to_exon(self) -> bool:
        return self.feature.kind == 101  # ord("e")

    @property
    def skip_makes_sense(self) -> bool:
        """If the SKIP in the segment matches some extremity of the feature and therefore can be interpreted as a splice event"""
        return (
            abs(self.feature.start - self.segment[0]) <= SPLIC_INACUR
            or abs(self.feature.end - self.segment[1]) <= SPLIC_INACUR
            if self.is_spliced
            else True
        )

    def __repr__(self) -> str:
        txt = "<SegmentMatch "
        if self.maps_to_intron:
            txt += "intron "
        if self.maps_to_exon:
            txt += "exon "
        if self.is_spliced:
            txt += "spliced"
        txt += f"\nSegmentPosition:{self.segment[0]}-{self.segment[1]} ({self.segment[1]-self.segment[0]+1}bp)"
        txt += f"\n{self.feature}\n>"
        return txt
