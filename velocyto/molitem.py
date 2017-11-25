from typing import *
import velocyto as vcy
import numpy as np


class Molitem:
    """ Keeps track of the alignments (to intervals of one gene) of all reads corresponding to the same molecule """

    __slots__ = ["gene", "has_some_spliced_read", "ivlhits"]

    def __init__(self, gene: vcy.Gene) -> None:
        self.gene = gene
        self.has_some_spliced_read = False  # type: bool
        self.ivlhits = np.zeros(gene.num_ivls(), dtype=bool)

    def mark_hit_ivls(self, matchivls: List[vcy.Interval], read_is_spliced: bool = False) -> None:
        """ Add info for the alignment(s) of one read """
        for matchivl in matchivls:
            self.ivlhits[matchivl.ivlidx] = True
        self.has_some_spliced_read = self.has_some_spliced_read or read_is_spliced

    def count(self, my_bcidx: int) -> None:
        """ Call after all reads have been processed to annotate this molecule on the gene """
        num_sure_exon_ivls = num_maybe_exon_ivls = 0
        # Loop trough the gene's intervals
        for hit, ivl in zip(self.ivlhits, self.gene.ivls):
            if hit:
                # if there is at least one evidence of sure intron mapping count as a unspliced molecule
                if ivl.is_sure_valid_intron:
                    self.gene.unspliced_mol_counts[my_bcidx] += 1
                    return
                num_sure_exon_ivls += ivl.is_sure_exon
                num_maybe_exon_ivls += ivl.is_maybe_exon
                # NOTE here we could keep track of ivl.is_sure_intron for the purpose of excluding them in the following condition
        if (num_sure_exon_ivls > 0) or self.has_some_spliced_read:  # NOTE: this condition might be reconsidered
            # alternative might be if (num_sure_exon_ivls > 0 and num_maybe_exon_ivls == 0) or self.has_some_spliced_read:  #
            # Either all reads hit to sure exons, or all reads hit to sure/maybe exon and some read is spliced
            self.gene.spliced_mol_counts[my_bcidx] += 1
        else:
            self.gene.ambiguous_mol_counts[my_bcidx] += 1
