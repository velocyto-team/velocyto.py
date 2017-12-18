from typing import *
import velocyto as vcy


class TranscriptModel:
    """A simple object representing a transcript model as a list of `vcy.Feature` objects
    """
    __slots__ = ["trid", "trname", "geneid", "genename", "chromstrand", "list_features"]
    
    def __init__(self, trid: str, trname: str, geneid: str, genename: str, chromstrand: str) -> None:
        self.trid = trid
        self.trname = trname
        self.geneid = geneid
        self.genename = genename
        self.chromstrand = chromstrand
        self.list_features: List[vcy.Feature] = []
    
    def __iter__(self) -> vcy.Feature:
        for i in self.list_features:
            yield i

    def __lt__(self, other: Any) -> bool:
        assert self.chromstrand == other.chromstrand, "`<`(.__lt__) not implemented for different chromosomes"
        return self.list_features[0].start < other.list_features[0].start

    def __gt__(self, other: Any) -> bool:
        assert self.chromstrand == other.chromstrand, "`>` (.__gt__) not implemented for different chromosomes"
        return self.list_features[0].start > other.list_features[0].start

    @property
    def start(self) -> int:
        """ NOTE: This should be accessed only after the creation of the transcript model is finished
        (i.e.) after append_exon has been called to add all the exons/introns
        """
        return self.list_features[0].start

    @property
    def end(self) -> int:
        """NOTE: This should be accessed only after the creation of the transcript model is finished
        (i.e.) after append_exon has been called to add all the exons/introns
        """
        return self.list_features[-1].end
        
    def ends_upstream_of(self, read: vcy.Read) -> bool:
        # one could consider to add TOLERANCE
        # note that ``self.list_features[-1]`` is the last exon if strand + and first exons for strand -
        return self.list_features[-1].end < read.pos

    def intersects(self, segment: Tuple[int, int], minimum_flanking: int=vcy.MIN_FLANK) -> bool:
        return (segment[-1] - minimum_flanking > self.start) and\
               (segment[0] + minimum_flanking < self.end)  # and ((segment[-1] - segment[0]) > minimum_flanking)

    def append_exon(self, exon_feature: vcy.Feature) -> None:
        """Append an exon and create an intron when needed

        Arguments
        ---------
        exon_feature: vcy.Feature
            A feature object represneting an exon to add to the transcript model.
        """
        exon_feature.transcript_model = self
        if len(self.list_features) == 0:
            # first/last exon
            self.list_features.append(exon_feature)
        else:
            # Some exon already exissted
            if self.chromstrand[-1] == "+":
                intron_number = self.list_features[-1].exin_no
            else:
                intron_number = self.list_features[-1].exin_no - 1
            self.list_features.append(vcy.Feature(start=self.list_features[-1].end + 1,
                                                  end=exon_feature.start - 1,
                                                  kind=ord("i"),
                                                  exin_no=intron_number,
                                                  transcript_model=self))
            self.list_features.append(exon_feature)

    def chop_if_long_intron(self, maxlen: int=vcy.LONGEST_INTRON_ALLOWED) -> None:
        """Modify a Transcript model choppin the 5' region upstram of a very long intron
        To avoid that extremelly long intron mask the counting of interal genes

        Arguments
        ---------
        maxlen: int, default=vcy.LONGEST_INTRON_ALLOWED
            transcript model tha contain one or more intronic interval of len == maxlen will be chopped

        Returns
        -------
        Nothing it will call `_remove_upstream_of` or `_remove_downstream_of` on the transcript model
        its name will be changed appending `_mod` to both trid and trname
         
        """
        long_feats = [i for i in self.list_features if len(i) > maxlen and i.kind == ord("i")]
        if len(long_feats):
            if self.chromstrand[-1] == "+":
                self._remove_upstream_of(long_feats[-1])
            else:  # self.chromstrand[-1] == "-"
                self._remove_downstream_of(long_feats[0])
            self.trid = self.trid + "_mod"
            self.trname = self.trname + "_mod"

    def _remove_upstream_of(self, longest_feat: vcy.Feature) -> None:
        tmp = []
        ec = 1
        ic = 1
        for feat in self.list_features:
            if feat > longest_feat:
                if feat.kind == ord("e"):
                    feat.exin_no = ec
                    ec += 1
                    tmp.append(feat)
                elif feat.kind == ord("i"):
                    feat.exin_no = ic
                    ic += 1
                    tmp.append(feat)
        self.list_features = tmp

    def _remove_downstream_of(self, longest_feat: vcy.Feature) -> None:
        tmp = []
        ec = 1
        ic = 1
        for feat in self.list_features[::-1]:
            if feat < longest_feat:
                if feat.kind == ord("e"):
                    feat.exin_no = ec
                    ec += 1
                    tmp.append(feat)
                elif feat.kind == ord("i"):
                    feat.exin_no = ic
                    ic += 1
                    tmp.append(feat)
        self.list_features = tmp[::-1]

    def __repr__(self) -> str:
        list_feats = '-'.join(f"{chr(i.kind)}{i.exin_no}" for i in self.list_features)
        return f"<TrMod {self.trid}\t{list_feats}\tat {hex(id(self))}>"
