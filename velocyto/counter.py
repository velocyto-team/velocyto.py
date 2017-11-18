import re
import sys
import gzip
import logging
from typing import *
from collections import defaultdict
import velocyto as vcy
import pysam
import numpy as np

tagpat = re.compile("NH:i:1.+CB:Z:([ACGT]+).+UB:Z:([ACGT]+)")


class ExInCounter:
    """ Main class to do the counting of introns and exons """
    def __init__(self, valid_bcs2idx: Dict[str, int]=None) -> None:
        """ valid_bcs2idx is a dict mapping valid cell barcodes to unique arbitrary numeric indexes used in the output arrays """
        if valid_bcs2idx is None:
            self.valid_bcs2idx: Dict[str, int] = {}
            self.filter_mode = False
            self.counter = 0
        else:
            self.valid_bcs2idx = valid_bcs2idx
            self.filter_mode = True
        self.mask_ivls_by_chromstrand = defaultdict(list)  # type: Dict[str, List]
        # NOTE: by using a default dict and not logging access to keys that do not exist, we might miss bugs!!!
        self.genes = []  # type: List
        self.ivls_by_chromstrand = defaultdict(list)  # type: Dict[str, List]

    @property
    def bclen(self) -> int:
        try:
            return len(next(iter(self.valid_bcs2idx.keys())))
        except StopIteration:
            return None

    def process(self, ivlfile: str, samfile: str) -> None:
        """ Performs all steps in the right order """
        self.read_genes(ivlfile)
        self.mark_up_introns(samfile)
        self.count(samfile)

    @staticmethod
    def parse_cigar_tuple(cigartuples: List[Tuple], pos: int) -> Tuple[List[Tuple[int, int]], bool, int, int]:
        segments = []
        spliced = False
        clip5 = clip3 = 0
        p = pos
        for operation_id, length in cigartuples:
            if operation_id == 0:  # vcy.CIGAR[operation_id] == "BAM_CMATCH"
                segments.append((p, p + length))
                p += length
            elif operation_id == 3:  # A splice || vcy.CIGAR[operation_id] == 'BAM_CREF_SKIP'
                spliced = True
                p += length
            elif operation_id == 2:  # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
                p += length
            elif operation_id == 4:  # bases at 5' or 3' are NOT part of the alignment || vcy.CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
                if p == pos:
                    clip5 = length  # At start of alignment
                else:
                    clip3 = length  # Must be at end of alignment vcy.CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
                p += length
            elif operation_id == 1:  # An insertion BAM_CINS
                pass  # maybe there should be a limit lenght accepted?
            elif operation_id == 5:  # BAM_CHARD_CLIP
                logging.warn("Hard clip was encountered! All mapping are assumed soft clipped")

        return segments, spliced, clip5, clip3

    def peek(self, samfile: str, lines: int=30) -> None:
        """Peeks into the samfile to determine if it is a cellranger or dropseq file
        """
        logging.debug(f"Peeking into {samfile}")
        fin = pysam.AlignmentFile(samfile)  # type: pysam.AlignmentFile
        cellranger: int = 0
        dropseq: int = 0
        failed: int = 0
        for i, read in enumerate(fin):
            if read.is_unmapped:
                continue
            if read.has_tag("CB") and read.has_tag("UB"):
                cellranger += 1
            elif read.has_tag("XC") and read.has_tag("XM"):
                dropseq += 1
            else:
                logging.warn(f"Not found cell and umi barcode in entry {i} of the bam file")
                failed += 1
            if cellranger > lines:
                self.cellbarcode_str = "CB"
                self.umibarcode_str = "UB"
                break
            elif dropseq > lines:
                self.cellbarcode_str = "XC"
                self.umibarcode_str = "XM"
                break
            elif failed > 5 * lines:
                raise IOError("The bam file does not contain cell and umi barcodes appropriatelly formatted")
            else:
                pass
        fin.close()

    def iter_unique_alignments(self, samfile: str, yield_line: bool=False) -> Iterable:
        """Iterates over the bam/sam file
        """
        self.peek(samfile, lines=10)
        logging.debug(f"Reading {samfile}")
        fin = pysam.AlignmentFile(samfile)  # type: pysam.AlignmentFile
        for i, read in enumerate(fin):
            if read.is_unmapped:
                continue
            try:
                bc = read.get_tag(self.cellbarcode_str).rstrip("-1")  # NOTE: this rstrip is relevant only for cellranger, should not cause trouble in Dorpseq
                umi = read.get_tag(self.umibarcode_str)
            except KeyError:
                continue  # NOTE: Here errors could go unnoticed
            if bc not in self.valid_bcs2idx:
                if self.filter_mode:
                    continue
                else:
                    self.valid_bcs2idx[bc] = self.counter
                    self.counter += 1
            strand = '-' if read.is_reverse else '+'
            chrom = fin.get_reference_name(read.rname)  # this is return a string otherwise read.rname for the integer
            if chrom.startswith('chr'):
                # I have to deal with incongruences with the .gft (what is cellranger doing???)
                # NOTE Why this happens?
                if "_" in chrom:
                    chrom = chrom.split("_")[1]
                else:
                    chrom = chrom[3:]
                    if chrom == "M":
                        chrom = "MT"
            pos = read.pos
            segments, spliced, clip5, clip3 = self.parse_cigar_tuple(read.cigartuples, pos)
            if segments == []:
                logging.debug("No segments in read:%s" % read.qname)
            if i % 10000000 == 0:
                logging.debug(f"Read first {i // 1000000} million mapped reads")
            if yield_line:
                yield vcy.Read(bc, umi, chrom, strand, pos, read.cigarstring, segments, clip5, clip3, spliced), read.tostring(fin)
            else:
                yield vcy.Read(bc, umi, chrom, strand, pos, read.cigarstring, segments, clip5, clip3, spliced)
        fin.close()

    def read_repeats(self, gtf_file: str) -> int:
        """ Read masked repeats """
        fin = open(gtf_file)
        for line in fin:
            fields = line.rstrip().split('\t')
            chrom, _, _, start, end, _, strand, _, _ = fields
            # Properly format the chromosome to match how it is parsed in the bam file
            if "_" in chrom:
                chrom = chrom.split("_")[1]
            else:
                chrom = chrom[3:]
                if chrom == "M":
                    chrom = "MT"
            ivlstart = int(start)
            ivlend = int(end)
            chromstrand = chrom + strand
            ivl = vcy.Interval(ivlstart, ivlend, None, 0, "mm10_rmsk", None)
            self.mask_ivls_by_chromstrand[chromstrand].append(ivl)

        # Make sure they are sorted (they should be becouse of the way the files where genrated)
        n = 0
        for chromstrand, ivls in self.mask_ivls_by_chromstrand.items():
            ivls.sort()  # relies on the __lt__ method of vcy.Interval
            n += len(ivls)
        return n

    def read_genes(self, ivlfile: str) -> int:
        """ Read genes and intervals """
        nbcs = len(self.valid_bcs2idx)
        for line in open(ivlfile):
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            trname, trid, genename, geneid, chrom, strand = fields[:6]
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            chromstrand = chrom + strand
            ivls = []  # type: List
            gene = vcy.Gene(genename, geneid, chrom, strand, nbcs)
            self.genes.append(gene)
            genestart, geneend = 2**30, 0
            for i in range(6, len(fields)):
                field = fields[i]
                span, ivltype, other_exons, other_introns = field.split(':')
                ivlstart, ivlend = span.split('-')
                ivlstart, ivlend = int(ivlstart), int(ivlend)
                genestart = min(genestart, ivlstart)
                geneend = max(geneend, ivlend)
                is_last_3exon_flag = (i == len(fields) - 1 and strand == "+") or (i == 6 and strand == "-")
                ivl = vcy.Interval(ivlstart, ivlend, gene, len(ivls), ivltype, is_last_3exon_flag)
                self.ivls_by_chromstrand[chromstrand].append(ivl)
                ivls.append(ivl)
            gene.set_range(genestart, geneend)
            gene.set_ivls(ivls)
        n = 0
        for chromstrand, ivls in self.ivls_by_chromstrand.items():
            ivls.sort()  # relies on the __lt__ method of vcy.Interval
            n += len(ivls)
        return n

    def mark_up_introns(self, samfile: str, intron_validation: str) -> None:
        """ Mark up introns that have reads across junctions to exons ?and count last exon reads? """
        # Read the file
        currchrom = ""
        for r in self.iter_unique_alignments(samfile):
            # Don't consider spliced reads (exonic) in this step
            if r.is_spliced():
                continue

            # When the chromosome mapping of the read changes, change interval index.
            if r.chrom != currchrom:
                logging.debug(f"Marking up chromosome {r.chrom}")
                currchrom = r.chrom
                if currchrom + '+' not in self.ivls_by_chromstrand:
                    logging.warn(f"The .bam file refers to a chromosome '{currchrom}' not present in the annotation (.gtf) file")
                iif = vcy.IntervalsIndex(self.ivls_by_chromstrand[currchrom + '+'])
                iir = vcy.IntervalsIndex(self.ivls_by_chromstrand[currchrom + '-'])
            
            # Consider the correct strand
            ii = iif if r.strand == '+' else iir

            # Look for overlap between the intervals and the read
            matchgenes, matchivls = ii.find_overlapping_ivls(r)
            if len(matchgenes) != 1:  # No match, multiple-matches, chimera, or intron read accidentally spanning another genes' exon zipped between the first genes' exons
                continue

            matchgene = matchgenes.pop()
            matchgene.add_read_stats(r)  # NOTE: this might be removed if the stats output won't be used... but not sure
            
            for matchivl, matchtype in matchivls.items():
                if matchtype == vcy.MATCH_INSIDE:
                    matchivl.add_read_inside()
                if matchtype & vcy.MATCH_OVER5END:  # NOTE: should this be elif? or it makes sense to allow both?
                    matchivl.add_read_spanning5end()
                if matchtype & vcy.MATCH_OVER3END:  # Added recently
                    matchivl.add_read_spanning3end()

        # Validate the intron intervals on the basis of occurrences of Exon/Intron spanning reads
        for gene in self.genes:
            gene.validate_intron_ivls(rule=intron_validation)  # NOTE: this condition might be stricter than we want

        # NOTE: did we properly deal with the reads added to the genes so that self.count() does not get confused?

    def count(self, samfile: str, sam_output: Tuple[Any, Any, Any, Any, Any]=None) -> None:
        """ Do the counting of molecules"""
        molitems = {}  # type: Dict[str, vcy.Molitem]  # keys are the cell_barcode + UMI, values are Molitem objects
        currchrom = ""

        if sam_output is not None:
            f_sure_introns, f_sure_exon, f_maybe_exon, f_others, f_chimeras = sam_output

        # Loop through the aligment of the samfile
        for r in self.iter_unique_alignments(samfile, yield_line=sam_output is not None):
            if sam_output is not None:
                r, line = r

            # When the chromosome mapping of the read changes, change interval index.
            if r.chrom != currchrom:
                currchrom = r.chrom
                logging.debug("Counting on chr %s" % currchrom)
                iif = vcy.IntervalsIndex(self.ivls_by_chromstrand[currchrom + '+'])
                iir = vcy.IntervalsIndex(self.ivls_by_chromstrand[currchrom + '-'])
                iimaskf = vcy.IntervalsIndex(self.mask_ivls_by_chromstrand[currchrom + '+'])
                iimaskr = vcy.IntervalsIndex(self.mask_ivls_by_chromstrand[currchrom + '-'])
                if not (currchrom + '+') in self.mask_ivls_by_chromstrand or not (currchrom + '-') in self.mask_ivls_by_chromstrand:
                    logging.debug(f"{currchrom} is not repeat masked.")

            # Consider the correct strand
            ii = iif if r.strand == '+' else iir
            iim = iimaskf if r.strand == '+' else iimaskr

            # Check if read is fully inside a repeat, in that case skip it
            if iim.has_ivls_enclosing(r):
                continue

            # Look for overlap between the intervals and the read
            matchgenes, matchivls = ii.find_overlapping_ivls(r)
            # logging.debug("matchgenes : %s" % matchgenes)
            if len(matchgenes) != 1:  # No match, multiple-matches, chimera, or intron read accidentally spanning another genes' exon zipped between the first genes' exons
                if sam_output is not None:
                    logging.debug("Chimera inside read or no matching gene")
                    f_chimeras.write(line)
                continue
            
            # Count molecules by dealing with barcodes
            bcumi = r.bc + r.umi
            matchgene = matchgenes.pop()  # NOTE: how does it compare to matchgenes[0]? Is there some gc issue?

            try:
                # Check if the read's gene is the same for the same molecule (bc+umi)
                molitem = molitems[bcumi]
                if molitem.gene != matchgene:
                    # In that case ignore it
                    if sam_output is not None:
                        logging.debug("Chimera between different reads of same molecule")
                        f_chimeras.write(line)
                    continue
            except KeyError:
                # The bc/umi combination was seen for the first time, create a Molecule Item
                molitem = vcy.Molitem(matchgene)
                molitems[bcumi] = molitem

            # Add info about which interval was hit and whether it gives evidence of splicing
            molitem.mark_hit_ivls(matchivls.keys(), r.is_spliced())
            
            if sam_output:
                # Resolve what the read is for the output sam files
                # NOTE: It would be nice to be able to output different sam files by molecule type other than read type
                cum_is_maybe_exon = []
                cum_is_sure_exon = []
                cum_is_sure_intron = []
                for k, v in matchivls.items():
                    cum_is_maybe_exon.append(k.is_maybe_exon)
                    cum_is_sure_exon.append(k.is_sure_exon)
                    cum_is_sure_intron.append(k.is_sure_intron)
                if np.any(cum_is_sure_intron):
                    f_sure_introns.write(line)
                elif np.any(cum_is_sure_exon):
                    f_sure_exon.write(line)
                elif np.any(cum_is_maybe_exon):
                    f_maybe_exon.write(line)

        # After the whole file has been read, do the actual counting
        for bcumi, molitem in molitems.items():
            bc = bcumi[:self.bclen]  # extract the bc part from the bc+umi
            bcidx = self.valid_bcs2idx[bc]
            molitem.count(bcidx)
