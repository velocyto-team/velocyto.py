from typing import *
import velocyto as vcy


class Transcript:
    """Object rappresenting a transcript"""
    __slots__ = ["trname", "trid", "genename", "geneid", "chrom", "strand", "_exons", "exons_sorted"]
    
    def __init__(self, trname: str, trid: str, genename: str, geneid: str, chrom: str, strand: str) -> None:
        """Object rappresenting a transcript
        
        Attributes
        ----------
        trname: str, trid: str, genename: str, geneid: str, chrom: str, strand: str

        Methods
        -------
        get_chromstrand, add_exon, fuse, sorted_exons, sorted_ivls
        get_start, get_end, get_chrom, get_strand, is_fw, get_trname, get_genename, get_trid
        spans_over
        
        Note
        ----
        It only stores the exons, the introns are inferred as exon_n_end->exon_n+1_start
        """
        self.trname = trname
        self.trid = trid
        self.genename = genename
        self.geneid = geneid
        self.chrom = chrom
        self.strand = strand
        self._exons = []  # type: List[List[int]]

    def __lt__(self, other: Any) -> bool:
        if self.get_chrom() != other.get_chrom():
            return self.get_chrom() < other.get_chrom()
        elif self.get_strand() != other.get_strand():
            return self.get_strand() < other.get_strand()
        elif self.get_start() != other.get_start():
            return self.get_start() < other.get_start()
        else:
            return self.get_end() < other.get_end()

    def get_chromstrand(self) -> str:
        return self.chrom + self.strand

    def add_exon(self, start: int, end: int, fuse: bool=False) -> None:
        self._exons.append([start, end])
        self.exons_sorted = False
        if fuse:
            self.fuse()
    
    def fuse(self) -> None:
        self._exons.sort()
        self.exons_sorted = True
        i = 1
        while i < len(self._exons):
            if self._exons[i][0] <= self._exons[i - 1][1]:
                self._exons[i - 1][1] = self._exons[i][1]
                del self._exons[i]
            else:
                i += 1

    def sorted_exons(self) -> List[Any]:
        self._exons.sort()
        self.exons_sorted = True
        return self._exons
    
    def sorted_ivls(self) -> Iterable:
        """Returns sorted intervals of both exons and introns using only the exons
        """
        if len(self._exons) > 0:
            exonstart, exonend = self._exons[0][0:2]
            yield [exonstart, exonend, vcy.EXON]
            for i in range(1, len(self._exons)):
                exonstart = self._exons[i][0]
                yield [exonend, exonstart, vcy.INTRON]
                exonend = self._exons[i][1]
                yield [exonstart, exonend, vcy.EXON]

    def spans_over(self, interval: Tuple[int, int]) -> bool:
        return (self.get_start() < interval[-1]) and (self.get_end() > interval[0])

    def get_start(self) -> int:
        return self.sorted_exons()[0][0]

    def get_end(self) -> int:
        return self.sorted_exons()[-1][1]

    def get_chrom(self) -> str:
        return self.chrom

    def get_strand(self) -> str:
        return self.strand

    def is_fw(self) -> bool:
        return self.strand == '+'

    def get_trname(self) -> str:
        return self.trname

    def get_genename(self) -> str:
        return self.genename

    def get_trid(self) -> str:
        return self.trid

    def get_geneid(self) -> str:
        return self.geneid

    def __str__(self) -> Any:
        return "Transcript(TrID=%s TrName=%s GeneID=%s GeneName=%s Chrom=%s Strand=%s\n\tExons=%s)\ni\tIvls=%s" % \
            (self.trid, self.trname, self.geneid, self.genename, self.chrom, self.strand,
             ",".join("%s-%s" % (s[0], s[1]) for s in self.sorted_exons()),
             ",".join("%s-%s:%s" % (s[0], s[1], s[2]) for s in self.sorted_ivls()))
