from numpy import zeros
import array
import logging
from typing import *
import velocyto as vcy


class Gene:

    __slots__ = ["genename", "geneid", "chrom", "strand", "lastexon_end_pos", "spliced_mol_counts",
                 "ambiguous_mol_counts", "unspliced_mol_counts", "deduced_tr_end", "start", "end",
                 "read_start_counts_from_locus_end", "ivls", "ivlinside_read_counts", "ivljunction5_read_counts", "ivljunction3_read_counts"]
                 
    def __init__(self, genename: str, geneid: str, chrom: str, strand: str, nbcs: int, genestart: int=0, geneend: int=0) -> None:
        self.genename = genename
        self.geneid = geneid
        self.chrom = chrom
        self.strand = strand
        self.lastexon_end_pos = 0  # Relative to annotated (+ means after)
        self.spliced_mol_counts = zeros(shape=nbcs, dtype=vcy.LOOM_NUMERIC_DTYPE)  # NOTE: this is uint16
        self.ambiguous_mol_counts = zeros(shape=nbcs, dtype=vcy.LOOM_NUMERIC_DTYPE)
        self.unspliced_mol_counts = zeros(shape=nbcs, dtype=vcy.LOOM_NUMERIC_DTYPE)
        self.deduced_tr_end = 0
        self.start = genestart
        self.end = geneend
        self.read_start_counts_from_locus_end = array.array('H')

    def __lt__(self, other: Any) -> bool:
        if self.get_chrom() != other.get_chrom():
            return self.get_chrom() < other.get_chrom()
        elif self.get_strand() != other.get_strand():
            return self.get_strand() < other.get_strand()
        elif self.get_start() != other.get_start():
            return self.get_start() < other.get_start()
        else:
            return self.get_end() < other.get_end()

    def get_start(self) -> int:
        return self.start

    def get_end(self) -> int:
        return self.end

    def get_chrom(self) -> str:
        return self.chrom

    def get_strand(self) -> str:
        return self.strand

    def get_spliced_mol_counts(gene, bcidx: int) -> int:
        """ Molecule counts that are exonic in all transcript models [ordered by barcode index]"""
        return gene.spliced_mol_counts[bcidx]

    def get_ambiguous_mol_counts(gene, bcidx: int) -> int:
        """ Molecule counts that are exonic in some but not all transcript models [ordered by barcode index]"""
        return gene.ambiguous_mol_counts[bcidx]

    def get_unspliced_mol_counts(gene, bcidx: int) -> int:
        """ Molecule counts that are intronic by all transcript models [ordered by barcode index]"""
        return gene.unspliced_mol_counts[bcidx]

    def set_range(self, start: int, end: int) -> None:
        self.start, self.end = start, end

    def set_ivls(self, ivls: List) -> None:
        self.ivls = ivls
        self.ivlinside_read_counts = zeros(shape=len(ivls), dtype="uint32")  # array.array('I', [0] * len(ivls))
        self.ivljunction5_read_counts = zeros(shape=len(ivls), dtype="uint32")  # array.array('I', [0] * len(ivls))
        self.ivljunction3_read_counts = zeros(shape=len(ivls), dtype="uint32")  # added recently

    def num_ivls(self) -> int:
        return len(self.ivls)

    def get_lastexon_counts(self) -> Tuple[int, int]:
        """ return the tuple: (read count spanning most 3' intron-exon junction, read count on most 3' exon) """
        result = 0, 0
        lastidx = len(self.ivls) - 1
        i, step, stop = (0, 1, lastidx) if self.strand == '-' else (lastidx, -1, 0)
        if self.ivls[i].is_sure_exon:
            lastexon_count = self.ivlinside_read_counts[i]
            while i != stop and self.ivls[i + step].is_sure_exon:
                lastexon_count += self.ivljunction5_read_counts[i]
                i += step
                lastexon_count += self.ivlinside_read_counts[i]
            if i != stop and self.ivls[i + step].is_sure_intron:
                lastjunction_count = self.ivljunction5_read_counts[i]
                result = lastjunction_count, lastexon_count
        return result

    def get_lastexon_length(self) -> int:
        i = len(self.ivls) - 1 if self.strand == '+' else 0
        return self.ivls[i].end - self.ivls[i].start

    def get_last_3prime_exon_interval(self) -> int:
        i = len(self.ivls) - 1 if self.strand == '+' else 0
        return self.ivls[i]

    def validate_intron_ivls(self, rule: str="strict") -> None:
        """ Annotate the introns that can safely be used for intron molecule counting =>
        either flanking exon is sure an exon and has >= 1 read spanning junction 
        
        Now supports two different heursitics

        "strict" and "permissive"
        """
        for i in range(1, len(self.ivls) - 2):
            if self.ivls[i].is_sure_intron:
                if rule == "strict":
                    self.ivls[i].is_sure_valid_intron = \
                        self.ivls[i - 1].is_sure_exon and self.ivljunction5_read_counts[i] > 0 or \
                        self.ivls[i + 1].is_sure_exon and self.ivljunction5_read_counts[i + 1] > 0
                elif rule == "permissive":
                    self.ivls[i].is_sure_valid_intron = True
                else:
                    raise IOError(f"rule {rule} is not supported")


    def add_read_stats(self, read: vcy.Read) -> None:
        if self.strand == '+':
            self.deduced_tr_end = max(self.deduced_tr_end, read.end() + read.clip3)
            dist_from_end = self.end - read.start()
        else:
            dist_from_end = read.end() - self.start
            if self.deduced_tr_end == 0:
                self.deduced_tr_end = read.start() - read.clip5
            else:
                self.deduced_tr_end = min(self.deduced_tr_end, read.start() - read.clip5)
        if dist_from_end >= len(self.read_start_counts_from_locus_end):
            self.read_start_counts_from_locus_end.extend([0] * (dist_from_end + 100))
        c = self.read_start_counts_from_locus_end[dist_from_end]
        if c < vcy.MAX_USHORT:
            self.read_start_counts_from_locus_end[dist_from_end] = c + 1

    def get_deduced_tr_end(self) -> int:
        return self.deduced_tr_end

    def get_tr_end(self) -> int:
        return self.start if self.strand == '-' else self.end
