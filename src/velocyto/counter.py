import contextlib
import gzip
import random
import re
import string
from collections import Counter, OrderedDict, defaultdict
from itertools import chain
from pathlib import Path
from typing import DefaultDict, Iterable

import h5py
import numpy as np
import pysam
import scipy as sp
from loguru import logger
from tqdm.auto import tqdm

from velocyto.constants import LONGEST_INTRON_ALLOWED, LOOM_NUMERIC_DTYPE, PATCH_INDELS, PLACEHOLDER_UMI_LEN
from velocyto.feature import Feature
from velocyto.gene_info import GeneInfo
from velocyto.indexes import FeatureIndex
from velocyto.logic import Logic
from velocyto.molitem import Molitem
from velocyto.read import Read
from velocyto.transcript_model import TranscriptModel

import better_exceptions

better_exceptions.hook()


class ExInCounter:
    """Main class to do the counting of introns and exons"""

    def __init__(
        self,
        sampleid: str,
        logic: Logic,
        valid_bcset: set[str] = None,
        umi_extension: str = "no",
        onefilepercell: bool = False,
        dump_option: str = "0",
        outputfolder: str = "./",
        loom_numeric_dtype: str = LOOM_NUMERIC_DTYPE,
    ) -> None:
        self.outputfolder = outputfolder
        self.sampleid = sampleid
        self.loom_numeric_dtype = loom_numeric_dtype
        self.logic = logic()
        # NOTE: maybe there shoulb be a self.logic.verify_inputs(args) step at the end of init
        if valid_bcset is None:
            self.valid_bcset: set[str] = set()
            self.filter_mode = False
            self.counter = 0
        else:
            self.valid_bcset = valid_bcset
            self.filter_mode = True
        self.annotations_by_chrm_strand: dict[str, dict[str, TranscriptModel]] = {}
        self.mask_ivls_by_chromstrand = defaultdict(list)  # type: dict[str, list]
        self.geneid2ix: dict[str, int] = {}
        self.genes: dict[str, GeneInfo] = {}
        if umi_extension.lower() == "no":
            self.umi_extract = self._no_extension
        elif umi_extension.lower() == "chr":
            self.umi_extract = self._extension_chr
        elif umi_extension.lower() in {"gene", "gx"}:
            self.umi_extract = self._extension_Gene
        elif umi_extension.endswith("bp"):
            self.umi_bp = int(umi_extension[:-2])
            self.umi_extract = self._extension_Nbp
        elif umi_extension.lower() == "without_umi":
            self.umi_extract = self._placeolder_umi
        else:
            raise ValueError(f"umi_extension {umi_extension} is not allowed. Use `no`, `Gene` or `[N]bp`")
        if onefilepercell:
            self.cell_barcode_get = self._bam_id_barcode
        else:
            self.cell_barcode_get = self._normal_cell_barcode_get
        if self.logic.stranded:
            if self.logic.accept_discordant:
                self.count_cell_batch = self._count_cell_batch_stranded_discordant
            else:
                self.count_cell_batch = self._count_cell_batch_stranded
        else:
            self.count_cell_batch = self._count_cell_batch_non_stranded
        # NOTE: by using a default dict and not logger access to keys that do not exist, we might miss bugs!!!
        self.test_flag = None
        if dump_option[0] == "p":
            self.kind_of_report = "p"
            self.every_n_report = int(dump_option[1:])
        else:
            self.kind_of_report = "h"
            self.every_n_report = int(dump_option)
        self.report_state = 0
        self.cellbarcode_str = "NULL_BC"  # This value should never be used this is just to initialize it and detect if there are bugs downstream
        self.umibarcode_str = "NULL_UB"  # This value should never be used this is just to initialize it and detect if there are bugs downstream

    # NOTE: not supported anymore because now we support variable length barcodes
    # @property
    # def bclen(self) -> int:
    #     try:
    #         return len(next(iter(self.valid_bcset)))
    #     except StopIteration:
    #         return None

    @staticmethod
    def parse_cigar_tuple(cigartuples: list[tuple], pos: int) -> tuple[list[tuple[int, int]], bool, int, int]:
        segments = []
        hole_to_remove = set()
        ref_skip = False
        clip5 = clip3 = 0
        p = pos
        for i, (operation_id, length) in enumerate(cigartuples):
            match operation_id:
                case 0:  # CIGAR[operation_id] == "BAM_CMATCH"
                    segments.append((p, p + length - 1))
                    p += length
                case 1:  # An insertion BAM_CINS
                    if length <= PATCH_INDELS:
                        with contextlib.suppress(IndexError):
                            if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                                hole_to_remove.add(len(segments) - 1)
                case 2:  # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
                    if length <= PATCH_INDELS:
                        with contextlib.suppress(IndexError):
                            if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                                hole_to_remove.add(len(segments) - 1)
                    p += length
                case 3:  # A splice || CIGAR[operation_id] == 'BAM_CREF_SKIP'
                    ref_skip = True
                    p += length
                case 4:  # bases at 5' or 3' are NOT part of the alignment || CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
                    if p == pos:
                        clip5 = length  # At start of alignment
                    else:
                        clip3 = (
                            length  # Must be at end of alignment CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
                        )
                    p += length
                case 5:  # BAM_CHARD_CLIP
                    logger.warning("Hard clip was encountered! All mapping are assumed soft clipped")

        # Merge segments separated by small insertions and deletions
        for a, b in enumerate(sorted(hole_to_remove)):  # NOTE maybe sorted is not required realy
            segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])

        return segments, ref_skip, clip5, clip3

    def peek(self, bamfile: str, lines: int = 1000) -> None:
        """Peeks into the samfile to determine if it is a cellranger or dropseq file"""
        logger.debug(f"Peeking into {bamfile}")
        fin = pysam.AlignmentFile(bamfile)  # type: pysam.AlignmentFile
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
                logger.warning(f"Not found cell and umi barcode in entry {i} of the bam file")
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
                raise IOError(
                    "The bam file does not contain cell and umi barcodes appropriatelly formatted. If you are runnin UMI-less data you should use the -U flag."
                )
        fin.close()

    def peek_umi_only(self, bamfile: str, lines: int = 30) -> None:
        """Peeks for umi into the samfile to determine if it is a cellranger or dropseq file"""
        logger.debug(f"Peeking into {bamfile}")
        fin = pysam.AlignmentFile(bamfile)  # type: pysam.AlignmentFile
        cellranger: int = 0
        dropseq: int = 0
        failed: int = 0
        for i, read in enumerate(fin):
            if read.is_unmapped:
                continue
            if read.has_tag("UB"):
                cellranger += 1
            elif read.has_tag("XM"):
                dropseq += 1
            else:
                logger.warning(f"Not found cell and umi barcode in entry {i} of the bam file")
                failed += 1
            if cellranger > lines:
                self.umibarcode_str = "UB"
                break
            elif dropseq > lines:
                self.umibarcode_str = "XM"
                break
            elif failed > 5 * lines:
                raise IOError(
                    "The bam file does not contain umi barcodes appropriatelly formatted. If you are runnin UMI-less data you should use the -U flag."
                )
        fin.close()

    def _no_extension(self, read: pysam.AlignedSegment) -> str:
        return read.get_tag(self.umibarcode_str)

    def _extension_Nbp(self, read: pysam.AlignedSegment) -> str:
        return read.get_tag(self.umibarcode_str) + read.query_alignment_sequence[: self.umi_bp]

    def _extension_Gene(self, read: pysam.AlignedSegment) -> str:
        try:
            return f"{read.get_tag(self.umibarcode_str)}_{read.get_tag('GX')}"
        except KeyError:
            return f"{read.get_tag(self.umibarcode_str)}_withoutGX"

    def _placeolder_umi(self, read: pysam.AlignedSegment) -> str:
        return "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(PLACEHOLDER_UMI_LEN))

    def _extension_chr(self, read: pysam.AlignedSegment) -> str:
        return f"{read.get_tag(self.umibarcode_str)}_{read.rname}:{read.reference_start // 10000000}"

    def _normal_cell_barcode_get(self, read: pysam.AlignedSegment) -> str:
        return read.get_tag(self.cellbarcode_str).split("-")[0]

    def _bam_id_barcode(self, read: pysam.AlignedSegment) -> str:
        return f"{self._current_bamfile}"

    def iter_alignments(self, bamfiles: tuple[str], unique: bool = True, yield_line: bool = False) -> Iterable:
        """Iterates over the bam/sam file and yield Read objects

        Arguments
        ---------
        bamfiles: tuple[str]
            path to the bam files
        unique: bool
            yield only unique alignments
        yield_line: bool
            whether to yield the raw sam line

        Returns
        -------
        yields Read for each valid line of the bam file
        or a tuple (Read, sam_line) if ``yield_line==True``
        NOTE: At the file change it yields a `None`
        """
        # No peeking here, this will happen a layer above, and only once on the not sorted file! before it was self.peek(samfile, lines=10)

        # bamfile_name_seen: set[str] = set()
        counter_skipped_no_barcode = 0
        if Counter(bamfiles).most_common(1)[0][1] != 1:
            logger.warning("The bamfiles names are not unique. The full path to them will be used as unique identifier")
            use_basename = False
        else:
            use_basename = True
        for bamfile in bamfiles:
            self._current_bamfile = Path(bamfile).name if use_basename else str(bamfile)
            logger.debug(f"Reading {bamfile}")
            self.peek(bamfile)
            fin = pysam.AlignmentFile(bamfile)  # type: pysam.AlignmentFile
            for i, read in enumerate(tqdm(fin, desc="Ingesting bamfile reads with barcodes")):
                if i % 10000000 == 0:
                    logger.debug(f"Read first {i // 1000000} million reads")
                if read.is_unmapped:
                    continue
                # If unique parameter is set to True, skip non-unique alignments
                if unique and read.get_tag("NH") != 1:
                    continue
                try:
                    # NOTE: this rstrip is relevant only for cellranger, should not cause trouble in Dropseq
                    bc = self.cell_barcode_get(read)
                    umi = self.umi_extract(read)
                except KeyError as err:
                    if read.has_tag(self.cellbarcode_str) and read.has_tag(self.umibarcode_str):
                        raise KeyError(
                            f"Some errors in parsing the cell barcode has occurred {self.cellbarcode_str}, {self.umibarcode_str}\n{read}"
                        ) from err
                    counter_skipped_no_barcode += 1
                    continue  # NOTE: Here errors could go unnoticed
                if bc not in self.valid_bcset:
                    if self.filter_mode:
                        continue
                    else:
                        self.valid_bcset.add(bc)
                strand = "-" if read.is_reverse else "+"
                chrom = fin.get_reference_name(
                    read.rname
                )  # this is return a string otherwise read.rname for the integer
                if chrom.startswith("chr"):
                    # I have to deal with incongruences with the .gft (what is cellranger doing???)
                    # NOTE Why this happens?
                    if "_" in chrom:
                        chrom = chrom.split("_")[1]
                    else:
                        chrom = chrom[3:]
                        if chrom == "M":
                            chrom = "MT"
                pos = (
                    read.reference_start + 1
                )  # reads in pysam are always 0-based, but 1-based is more convenient to work with in bioinformatics
                segments, ref_skipped, clip5, clip3 = self.parse_cigar_tuple(read.cigartuples, pos)
                if segments == []:
                    logger.debug(f"No segments in read: {read.qname}")

                read_object = Read(bc, umi, chrom, strand, pos, segments, clip5, clip3, ref_skipped)
                if read_object.span > 3000000:  # Longest locus existing
                    logger.warning(f"Trashing read, too long span\n{read.tostring(fin)}")
                elif yield_line:
                    yield read_object, read.tostring(fin)
                else:
                    yield read_object
            fin.close()
            # NOTE Yielding None counts as a flag that the file read has been changed
            if yield_line:
                yield None, None
            else:
                yield None
        logger.debug(
            f"{counter_skipped_no_barcode} reads were skipped because no apropiate cell or umi barcode was found"
        )

    def read_repeats(self, gtf_file: str, tolerance: int = 5) -> None:
        """Read repeats and merge close ones into highly repetitive areas

        Arguments
        ---------
        gtf_file: str
            file to read
        tolerance: int, default=5
            if two repeats intervals to be masked are found closer than tolerance bases from each other they are fused in one bigger masked interval.
            Notice that in the downstream analysis only reads that are fall inside mask intervals are discarded

        Returns
        -------
        None
            why return anything?

        """
        # mask_ivls_by_chromstrand: dict[str, list[Feature]]
        #     A dictionary key: chromosome+strand value: list of features (repeat intervals)
        #     (The reference is returned but an internal attribure self.self.masked_by_chrm_strand is kept)
        # Example code to sort the gtf file
        # sorted_filename = gtffile.split(".")[0] + "_sorted.gtf"
        # logger.debug(f"Sorting by `sort -k1,1 -k7,7 -k4,4n {gtffile} > {sorted_filename}`")
        # with open(sorted_filename, "w") as f:
        #     p1 = subprocess.run(["sort", "-k1,1", "-k7,7", "-k4,4n", gtffile],
        #                         stdout=f)

        # Parse arguments
        logger.debug(f"Reading {gtf_file}, the file will be sorted in memory")

        # Read up skipping headers up to the first valid entry
        repeat_ivls_list: list[Feature] = []

        # fin = open(gtf_file)
        if gtf_file.suffix == ".gz":
            with gzip.open(gtf_file, "rb") as gtf:
                gtf_lines = [line.decode() for line in gtf if not line.decode().startswith("#")]
        else:
            with open(gtf_file, "r") as gtf:
                gtf_lines = [line for line in gtf if not line.startswith("#")]

        def sorting_key(entry: str) -> tuple[str, bool, int, str]:
            """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
            x = entry.split("\t")
            return (
                x[0],
                x[6] == "+",
                int(x[3]),
                entry,
            )  # The last element of the touple corresponds to the `last resort comparison`

        gtf_lines = sorted(gtf_lines, key=sorting_key)
        ###

        line = gtf_lines.pop(0)
        fields = line.rstrip().split("\t")
        (
            chrom,
            feature_class,
            feature_type,
            start_str,
            end_str,
            junk,
            strand,
            junk,
            tags,
        ) = fields
        # Removing chr from the chromosome name to uniform different formats of gtf files, taht might or might not have the prefix "chr"
        if chrom[:3].lower() == "chr":
            chrom = chrom[3:]
        start = int(start_str)
        end = int(end_str)
        chromstrand: str = chrom + strand

        # Set this tu the current entry
        curr_chrom: str = chrom
        # curr_feature_class: str = feature_class
        # curr_feature_type: str = feature_type
        curr_start: int = start
        curr_end: int = end
        curr_n: int = 1
        curr_strand: str = strand
        curr_tags: str = tags
        curr_chromstrand: str = chromstrand

        for line in tqdm(gtf_lines, desc="parsing GTF"):
            fields = line.rstrip().split("\t")
            (
                chrom,
                feature_class,
                feature_type,
                start_str,
                end_str,
                junk,
                strand,
                junk,
                tags,
            ) = fields
            # Removing chr from the chromosome name to uniform different formats of gtf files, that might or might not have the prefix "chr"
            if chrom[:3].lower() == "chr":
                chrom = chrom[3:]
            start = int(start_str)
            end = int(end_str)
            chromstrand = chrom + strand

            # On chromosome or strand change: save and clean after yourself
            if chromstrand != curr_chromstrand:
                self.mask_ivls_by_chromstrand[curr_chromstrand] = repeat_ivls_list
                repeat_ivls_list = []
                curr_chrom = chrom
                curr_strand = strand
                curr_chromstrand = curr_chrom + curr_strand

            # If there is no overlap/contiguity within a certain tolerance with the next
            if start > curr_end + tolerance:
                # Close and store the previus interval
                repeat_ivls_list.append(Feature(start=curr_start, end=curr_end, kind=ord("r"), exin_no=curr_n))
                # NOTE exin_no is being used to store curr_n, and kind could be used to store tag for some extra memory cost

                # Remember the newly parsed interval
                curr_start = start
                curr_end = end
                curr_n = 1  # number of original repeat regions included

                # this is extra information that right now is not stored in Feature (not to waste memory) since is not used
                # curr_feature_class = feature_class
                # curr_feature_type = feature_type
                curr_tags = tags
            else:
                # Extend the masked interval!
                curr_end = end
                curr_n += 1
                # this is extra information that right now is not saved (not to waste memory) since is not used
                gap = start - curr_end
                curr_tags = f"{curr_tags} gap {gap}; {tags}" if gap > 0 else curr_tags + tags

        n = 0
        for _, feature_list in self.mask_ivls_by_chromstrand.items():
            feature_list.sort()  # relies on the __lt__ method of Feature
            n += len(feature_list)

        logger.debug(f"Processed masked annotation .gtf and generated {n} intervals to mask!")

    def assign_indexes_to_genes(self, features: dict[str, TranscriptModel]) -> None:
        """Assign to each newly encoutered gene an unique index corresponding to the output matrix column ix"""
        logger.debug("Assigning indexes to genes")
        for trmodel in features.values():
            if trmodel.geneid in self.geneid2ix:
                self.genes[trmodel.geneid].start = min(self.genes[trmodel.geneid].start, trmodel.start)
                self.genes[trmodel.geneid].end = max(self.genes[trmodel.geneid].end, trmodel.end)
            else:
                self.geneid2ix[trmodel.geneid] = len(self.geneid2ix)
                self.genes[trmodel.geneid] = GeneInfo(
                    trmodel.genename,
                    trmodel.geneid,
                    trmodel.chromstrand,
                    trmodel.start,
                    trmodel.end,
                )

    def read_transcriptmodels(self, gtf_file: str) -> dict[str, dict[str, TranscriptModel]]:
        """Reads transcript models from a sorted .gtf file

        Arguments
        ---------
        gtf_file: str
            Path to the sorted gtf file

        Returns
        -------
        annotations_by_chrm_strand: dict[str, list[TrancriptModel]]
            A dictionary key: chromosome+strand value: list of trascript models
            (The reference is returned but an internal attribure self.annotations_by_chrm_strand is kept)

        There will exist an object Features for the same exon appearing in a different TranscriptModel. (his is desired)
        """

        # Define the regexes for the parsing
        regex_trid = re.compile(r'transcript_id "([^"]+)"')
        regex_trname = re.compile(r'transcript_name "([^"]+)"')
        regex_geneid = re.compile(r'gene_id "([^"]+)"')
        regex_genename = re.compile(r'gene_name "([^"]+)"')
        regex_exonno = re.compile(r'exon_number "*?([\w]+)')  # re.compile('exon_number "([^"]+)"')

        # Initialize containers
        # headerlines: list[str] = []

        if gtf_file.suffix == ".gz":
            with gzip.open(gtf_file, "rb") as gtf:
                gtf_lines = [line.decode() for line in gtf if not line.decode().startswith("#")]
        else:
            with open(gtf_file, "r") as gtf:
                gtf_lines = [line for line in gtf if not line.startswith("#")]

        def sorting_key(entry: str) -> tuple[str, bool, int, str]:
            """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
            x = entry.split("\t")
            return (
                x[0],
                x[6] == "+",
                int(x[3]),
                entry,
            )  # The last element of the touple corresponds to the `last resort comparison`

        gtf_lines = self.peek_and_correct(gtf_lines)
        gtf_lines = sorted(gtf_lines, key=sorting_key)

        curr_chromstrand = None
        features: dict[str, TranscriptModel] = OrderedDict()
        # Loop throug gtf file (assumes it is ordered)
        for nth_line, line in enumerate(gtf_lines):
            # Deal with headers
            # if line.startswith('#'):
            #     headerlines.append(line)
            #     continue

            fields = line.rstrip().split("\t")
            (
                chrom,
                feature_class,
                feature_type,
                start_str,
                end_str,
                junk,
                strand,
                junk,
                tags,
            ) = fields

            # Deal with possible incongruences between .gtf and bam file chromosome naming
            # We keep the name of the chromosome removing the chr from the name and removing part after the . if present
            if "chr" in chrom[:4]:
                chrom = chrom[3:]  # NOTE before it was chrom[3:].split(".")[0]
            # A new chromosome/strand encountered
            if chrom + strand != curr_chromstrand:
                if curr_chromstrand is not None:  # Every time with exception with first and the last chromosome
                    if chrom + strand in self.annotations_by_chrm_strand:
                        # NOTE this is not enough as a check but it will detect with few checks if file is not sorted at all
                        raise IOError(
                            "Genome annotation gtf file is not sorted correctly! Run the following command:\nsort -k1,1 -k7,7 -k4,4n -o [GTF_OUTFILE] [GTF_INFILE]"
                        )
                    else:
                        logger.debug(f"Done with {curr_chromstrand} [line {nth_line-1}]")
                    self.assign_indexes_to_genes(features)
                    self.annotations_by_chrm_strand[curr_chromstrand] = features
                    logger.debug(f"Seen {len(self.geneid2ix)} genes until now")
                features = OrderedDict()
                logger.debug(f"Parsing Chromosome {chrom} strand {strand} [line {nth_line}]")
                curr_chromstrand = chrom + strand

            if feature_type in ("exon"):
                trid = regex_trid.search(tags)[1]
                _trname_search = regex_trname.search(tags)
                trname = trid if _trname_search is None else _trname_search[1]
                geneid = regex_geneid.search(tags)[1]
                _genename_search = regex_genename.search(tags)
                genename = geneid if _genename_search is None else _genename_search[1]
                try:
                    exonno = regex_exonno.search(tags)[1]
                except AttributeError as err:
                    # NOTE: Don't try to release this constraint, velocyto relies on it for safe calculations! Rather make a utility script that does putative annotation separatelly.
                    raise IOError(
                        "The genome annotation .gtf file provided does not contain exon_number. `exon_number` is described as a mandatory field by GENCODE gtf file specification and we rely on it for easier processing"
                    ) from err
                start = int(start_str)
                end = int(end_str)
                chromstrand = chrom + strand

                try:
                    features[trid].append_exon(Feature(start=start, end=end, kind=ord("e"), exin_no=exonno))
                except KeyError:
                    features[trid] = TranscriptModel(
                        trid=trid,
                        trname=trname,
                        geneid=geneid,
                        genename=genename,
                        chromstrand=chromstrand,
                    )
                    features[trid].append_exon(Feature(start=start, end=end, kind=ord("e"), exin_no=exonno))

        # Do it for the last chromosome
        self.assign_indexes_to_genes(features)
        self.annotations_by_chrm_strand[curr_chromstrand] = features
        logger.debug(f"Done with {curr_chromstrand} [line {nth_line-1}]")

        logger.debug(
            f"Fixing corner cases of transcript models containg intron longer than {LONGEST_INTRON_ALLOWED//1000}Kbp"
        )
        # Fix corner cases of extremelly long introns ~1Mbp that would be masking genes that are found internally
        for tmodels_orddict in self.annotations_by_chrm_strand.values():
            for tm in tmodels_orddict.values():
                tm.chop_if_long_intron()  # Change it in place

        # Respect the sorting it had before
        # NOTE: not sure it is needed downstream anymore also not sure it guarantees exactly the same order
        for chromstrand in self.annotations_by_chrm_strand.keys():
            tmp = OrderedDict((i.trid, i) for i in sorted(self.annotations_by_chrm_strand[chromstrand].values()))
            self.annotations_by_chrm_strand[chromstrand] = tmp

        return self.annotations_by_chrm_strand

    def peek_and_correct(self, gtf_lines: list[str]) -> list[str]:
        """Look at the first 20 instances of a list of lines of a gtf file to dermine if exon number is specified as it should.
        If econ number is not contained it will infer the exon number sorting the list by lexicographic ordering tr_id, start, end

        Arguments
        ---------
        gtf_lines:
            a list of the lines of a gtf file

        Returns
        -------
        gtf_lines:
            the same list or the list corrected with added a exon number (filtered to contain only exons)
        """
        regex_exonno = re.compile(r'exon_number "*?([\w]+)')
        flag = False
        # Scan the first 500 lines for an occurrence of exon
        for lin in gtf_lines[:500]:
            (
                chrom,
                feature_class,
                feature_type,
                start_str,
                end_str,
                junk,
                strand,
                junk,
                tags,
            ) = lin.split("\t")
            if feature_type == "exon":
                exonno = regex_exonno.search(tags)
                if exonno is None:
                    flag = True

        if not flag:
            return gtf_lines
        logger.warning("The entry exon_number was not present in the gtf file. It will be infferred from the position.")
        regex_trid = re.compile('transcript_id "([^"]+)"')
        min_info_lines_minus: list[list] = []
        min_info_lines_plus: list[list] = []
        for lin in gtf_lines:
            (
                chrom,
                feature_class,
                feature_type,
                start_str,
                end_str,
                junk,
                strand,
                junk,
                tags,
            ) = lin.split("\t")
            if feature_type == "exon":
                try:
                    trid = regex_trid.search(tags)[1]
                except AttributeError as err:
                    raise AttributeError(f"transcript_id entry not found in line: {lin}") from err
                if strand == "-":
                    min_info_lines_minus.append([trid, int(start_str), int(end_str), lin])
                else:
                    min_info_lines_plus.append([trid, int(start_str), int(end_str), lin])

        min_info_lines_minus = sorted(min_info_lines_minus)
        min_info_lines_plus = sorted(min_info_lines_plus)
        current_trid = "None"
        exon_n = 1
        modified_lines_plus: list[str] = []
        for i in min_info_lines_plus:
            if current_trid != i[0]:
                current_trid = i[0]
                exon_n = 1
            else:
                exon_n += 1
            modified_lines_plus.append(f'{i[3][:-1]} exon_number "{exon_n}";\n')
        exon_n = 1
        modified_lines_minus: list[str] = []
        for i in min_info_lines_minus[::-1]:
            if current_trid != i[0]:
                current_trid = i[0]
                exon_n = 1
            else:
                exon_n += 1
            modified_lines_plus.append(f'{i[3][:-1]} exon_number "{exon_n}";\n')

        return modified_lines_plus + modified_lines_minus

    def mark_up_introns(self, bamfile: tuple[str], multimap: bool) -> None:
        """Mark up introns that have reads across exon-intron junctions

        Arguments
        ---------
        bamfile: tuple[str]
            path to the bam files to markup
        logic: Logic
            The logic object to use, changes in different techniques / levels of strictness
            NOTE: Right now it is not used

        Returns
        -------
        Nothing it just add to validation to the Features

        Note
        ----
        Situations not considered:
        # If an the exon is so short that is possible to get both exonA-exonB junction and exonB-intronB boundary in the same read
        """

        if not self.logic.perform_validation_markup:
            return
        # Since I support multiple files (Smart seq2 it makes sense here to load the feature indexes into memory)
        # NOTE this means that maybe I could do this once at a level above
        # NOTE if this is not done in count then I need to bring it before the if/else statement
        self.feature_indexes: DefaultDict[str, FeatureIndex] = defaultdict(FeatureIndex)
        for (
            chromstrand_key,
            annotions_ordered_dict,
        ) in tqdm(self.annotations_by_chrm_strand.items(), desc="Reading feature indices"):
            self.feature_indexes[chromstrand_key] = FeatureIndex(
                sorted(chain.from_iterable(annotions_ordered_dict.values()))
            )

        # VERBOSE: dump_list = []
        # Read the file
        currchrom = ""
        set_chromosomes_seen: set[str] = set()
        # TODO: feel like self.iter_alignments should not be a generator as we are currently running through it twice?
        for r in tqdm(self.iter_alignments(bamfile, unique=not multimap), desc="Ingesting reads"):
            # Don't consider spliced reads (exonic) in this step
            # NOTE Can the exon be so short that we get splicing and exon-intron boundary
            if r is None:
                logger.debug("r is None")
                # This happens only when there is a change of file
                currchrom = ""
                set_chromosomes_seen = set()
                # I need to reset indexes in position before the next file is restarted
                # NOTE this is far from optimal for lots of cells
                logger.debug("End of file. Reset index: start scanning from initial position.")
                for (
                    chromstrand_key,
                    _,
                ) in self.annotations_by_chrm_strand.items():
                    self.feature_indexes[chromstrand_key].reset()
                continue
            if r.is_spliced:
                # logger.debug(f"{r.is_spliced=}")
                continue

                # When the chromosome mapping of the read changes, change interval index.
            if r.chrom != currchrom:
                logger.debug(f"{r.chrom=}, {currchrom=}")
                if r.chrom in set_chromosomes_seen:
                    msg = f"Input .bam file should be chromosome-sorted. (Hint: use `samtools sort {bamfile}`)"
                    raise IOError(msg)
                set_chromosomes_seen.add(r.chrom)
                logger.debug(f"Marking up chromosome {r.chrom}")
                currchrom = r.chrom
                if f"{currchrom}+" not in self.annotations_by_chrm_strand:
                    logger.warning(
                        f"The .bam file refers to a chromosome '{currchrom}+' not present in the annotation (.gtf) file"
                    )
                    iif = FeatureIndex([])
                else:
                    iif = self.feature_indexes[f"{currchrom}+"]
                    logger.debug(f"{iif=}")
                if f"{currchrom}-" not in self.annotations_by_chrm_strand:
                    logger.warning(
                        f"The .bam file refers to a chromosome '{currchrom}-' not present in the annotation (.gtf) file"
                    )
                    iir = FeatureIndex([])
                else:
                    iir = self.feature_indexes[f"{currchrom}-"]
                    logger.debug(f"{iif=}")

            # Consider the correct strand
            ii = iif if r.strand == "+" else iir

            # VERBOSE: # Look for overlap between the intervals and the read
            # VERBOSE: dump_list += ii.mark_overlapping_ivls(r)
            ii.mark_overlapping_ivls(r)

            # VERBOSE: import pickle
            # VERBOSE: pickle.dump(dump_list, open("dump_mark_overlapping_ivls.pickle", "wb"))

    def count(
        self,
        bamfile: tuple[str],
        multimap: bool,
        cell_batch_size: int = 100,
        molecules_report: bool = False,
    ) -> tuple[dict[str, list[np.ndarray]], list[str]]:
        """Do the counting of molecules

        Arguments
        ---------
        bamfile: str
            path to the bam files to markup
        cell_batch_size: int, default = 50
            it defines whether to require or not exon-intron spanning read to consider an intron valid.

        Returns
        -------
        dict_list_arrays, cell_bcs_order

        Note
        ----
        The memory footprint could be reduced allocating scipy sparse arrays

        """

        # Initialize variables and containers
        self.cell_batch: set[str] = set()
        self.reads_to_count: list[Read] = []
        # self.cells_since_last_count = 0
        # Analysis is cell wise so the Feature Index swapping is happening often and it is worth to preload everything in memory
        # NOTE: for the features this was already done at markup time, maybe I should just reset them
        self.feature_indexes: DefaultDict[str, FeatureIndex] = defaultdict(FeatureIndex)
        for chromstrand_key, annotions_ordered_dict in tqdm(
            self.annotations_by_chrm_strand.items(), desc="Count molecules - feature indexes"
        ):
            self.feature_indexes[chromstrand_key] = FeatureIndex(
                sorted(chain.from_iterable(annotions_ordered_dict.values()))
            )

        self.mask_indexes: DefaultDict[str, FeatureIndex] = defaultdict(FeatureIndex)
        for chromstrand_key, annotions_list in tqdm(
            self.mask_ivls_by_chromstrand.items(), desc="Count molecules - mask indexes"
        ):
            self.mask_indexes[chromstrand_key] = FeatureIndex(annotions_list)  # This suould be sorted

        logger.debug(f"Features available for chromosomes : {list(self.feature_indexes)}")
        logger.debug(f"Mask available for chromosomes : {list(self.mask_indexes)}")

        # Before counting, report how many features where validated
        logger.debug("Summarizing the results of intron validation.")
        n_is_intron = 0
        n_is_intron_valid = 0
        unique_valid = set()
        for _, feature_index in tqdm(self.feature_indexes.items(), desc="Count molecules - validate introns"):
            for ivl in feature_index.ivls:
                if ivl.kind == ord("i"):
                    n_is_intron += 1
                if ivl.is_validated:
                    n_is_intron_valid += 1
                    unique_valid.add((ivl.start, ivl.end))
        logger.debug(
            f"Validated {n_is_intron_valid} introns (of which unique intervals {len(unique_valid)}) out of {n_is_intron} total possible introns (considering each possible transcript models)."
        )

        cell_bcs_order: list[str] = []
        dict_list_arrays: dict[str, list[np.ndarray]] = {layer_name: [] for layer_name in self.logic.layers}
        nth = 0
        # Loop through the aligment of the bamfile
        for r in tqdm(self.iter_alignments(bamfile, unique=not multimap), desc="Count molecules: count"):
            # if nth < 50:
            if (r is None) or (len(self.cell_batch) == cell_batch_size and r.bc not in self.cell_batch):
                # Perfrom the molecule counting
                nth += 1
                logger.debug(
                    f"Counting for batch {nth}, containing {len(self.cell_batch)} cells and {len(self.reads_to_count)} reads"
                )
                dict_layer_columns, list_bcs = self.count_cell_batch()

                # This is to avoid crazy big matrix output if the barcode selection is not chosen
                if not self.filter_mode:
                    logger.warning("The barcode selection mode is off, no cell events will be identified by <80 counts")
                    tot_mol = dict_layer_columns["spliced"].sum(0) + dict_layer_columns["unspliced"].sum(0)
                    cell_bcs_order += list(np.array(list_bcs)[tot_mol > 80])
                    for layer_name, layer_columns in dict_layer_columns.items():
                        dict_list_arrays[layer_name].append(layer_columns[:, tot_mol > 80])
                    logger.warning(f"{np.sum(tot_mol < 80)} of the barcodes where without cell")
                else:
                    # The normal case
                    cell_bcs_order += list_bcs
                    for layer_name, layer_columns in dict_layer_columns.items():
                        dict_list_arrays[layer_name].append(layer_columns)

                self.cell_batch = set()
                # Drop the counted reads (If there are no other reference to it) and reset the indexes to 0
                self.reads_to_count = []
                for (
                    chromstrand_key,
                    _,
                ) in self.annotations_by_chrm_strand.items():
                    self.feature_indexes[chromstrand_key].reset()
                for (
                    chromstrand_key,
                    _,
                ) in self.mask_ivls_by_chromstrand.items():
                    self.mask_indexes[chromstrand_key].reset()

            if r is not None:
                self.cell_batch.add(r.bc)
                self.reads_to_count.append(r)
        # NOTE: Since iter_allignments is yielding None at each file change (even when only one bamfile) I do not need the following
        # logger.debug(f"Counting molecule for last batch of {len(self.cell_batch)}, total reads {len(self.reads_to_count)}")
        # spliced, unspliced, ambiguous, list_bcs = self.count_cell_batch()
        # cell_bcs_order += list_bcs
        # list_spliced_arrays.append(spliced)
        # list_unspliced_arrays.append(unspliced)
        # list_ambiguous_arrays.append(ambiguous)
        # self.cell_batch = set()
        # self.reads_to_count = []
        logger.debug("Counting done!")
        return dict_list_arrays, cell_bcs_order

    def _count_cell_batch_stranded(self) -> tuple[dict[str, np.ndarray], list[str]]:
        """It performs molecule counting for the current batch of cells in the case of stranded method

        Returns
        -------
        dict_layers_columns: dict[str, np.ndarray]
            name_layer->np.ndarray of the batch
        idx2bc: list[str]
            list of barcodes

        NOTE This duplications of method is bad for code mantainance
        """
        molitems: DefaultDict[str, Molitem] = defaultdict(Molitem)
        # Sort similarly to what the sort linux command would do. (implemented using Read.__lt__)
        # NOTE NOTE!!!! Here I changed the way to sort because it was using the strand causing to skip a lot of reads in SmartSeq2
        self.reads_to_count.sort()
        # NOTE: I could start by sorting the reads by chromosome, strand, position but for now let's see if it is fast without doing do

        repeats_reads_count = 0
        for r in self.reads_to_count:
            # Consider the correct strand
            ii = self.feature_indexes[r.chrom + r.strand]
            iim = self.mask_indexes[r.chrom + r.strand]

            # Check if read is fully inside a masked region, in that case skip it
            if iim.has_ivls_enclosing(r):
                repeats_reads_count += 1  # VERBOSE
                continue

            # Look for overlap between the intervals and the read
            mappings_record = ii.find_overlapping_ivls(r)
            if len(mappings_record):
                # logger.debug("IN: Non empty mapping record")
                bcumi = f"{r.bc}${r.umi}"
                molitems[bcumi].add_mappings_record(mappings_record)
                # if len(molitems[bcumi].mappings_record):
                #     logger.debug("OUT: Non empty")
                # else:
                #     logger.debug("OUT: Empty")
        if repeats_reads_count > 0:
            logger.debug(
                f"{repeats_reads_count} reads not considered because fully enclosed in repeat masked regions"
            )  # VERBOSE
            # initialize np.ndarray
        shape = (len(self.geneid2ix), len(self.cell_batch))

        # dict_layers_columns: dict[str, np.ndarray] = {}
        dict_layers_columns: dict[str, sp.sparse.lil_array] = {
            layer_name: sp.sparse.lil_array(shape, dtype=self.loom_numeric_dtype) for layer_name in self.logic.layers
        }
        # for layer_name in self.logic.layers:
        #     # dict_layers_columns[layer_name] = np.zeros(shape, dtype=self.loom_numeric_dtype, order="C")
        #     dict_layers_columns[layer_name] = sp.sparse.lil_array(shape, dtype=self.loom_numeric_dtype)

        bc2idx: dict[str, int] = dict(zip(self.cell_batch, range(len(self.cell_batch))))
        # After the whole file has been read, do the actual counting
        failures = 0
        counter: Counter = Counter()
        for bcumi, molitem in molitems.items():
            bc = bcumi.split("$")[0]  # extract the bc part from the bc+umi
            bcidx = bc2idx[bc]
            if rcode := self.logic.count(molitem, bcidx, dict_layers_columns, self.geneid2ix):
                failures += 1
                counter[rcode] += 1
                # before it was molitem.count(bcidx, spliced, unspliced, ambiguous, self.geneid2ix)
        if failures > (0.25 * len(molitems)):
            logger.warning(f"More than 20% ({(100*failures / len(molitems)):.1f}%) of molitems trashed, of those:")
            if (x := (100 * counter[1] / len(molitems))) > 0.0:
                logger.warning(f"A situation where many genes were compatible with the observation in {x:.1f} cases")
            if (y := (100 * counter[2] / len(molitems))) > 0.0:
                logger.warning(f"No gene is compatible with the observation in {y:.1f} cases")
            if (z := (100 * counter[3] / len(molitems))) > 0.0:
                logger.warning(f"Observation compatible with more genes {z:.1f} of the cases")
            if (w := (100 * counter[4] / len(molitems))) > 0.0:
                logger.warning(f"Situation that were not described by the logic in the {w:.1f} of the cases")

        if self.every_n_report and ((self.report_state % self.every_n_report) == 0):
            if self.kind_of_report == "p":
                import pickle

                first_cell_batch = next(iter(molitems.keys())).split("$")[0]
                if not (x := Path("pickle_dump")).exists():
                    x.mkdir()
                pickle.dump(
                    molitems,
                    open(f"pickle_dump/molitems_dump_{first_cell_batch}.pickle", "wb"),
                )
                pickle.dump(
                    self.reads_to_count,
                    open(f"pickle_dump/reads_to_count{first_cell_batch}.pickle", "wb"),
                )
            else:
                if not (x := Path(self.outputfolder).joinpath("dump")).exists():
                    x.mkdir()
                f = h5py.File(
                    Path(self.outputfolder).joinpath("dump", f"{self.sampleid}.hdf5")
                )  # From the docs: Read/write if exists, create otherwise (default)

                if "info/tr_id" not in f:
                    logger.warning(
                        "The hdf5 report is less accurate in reporting exactly all the information than the pickle."
                    )
                    info_tr_id = []
                    info_features_gene = []
                    info_is_last3prime = []
                    info_is_intron = []
                    info_start_end = []
                    info_exino = []
                    info_strandplus = []
                    info_chrm = []
                    for _k, v_dict_tm in self.annotations_by_chrm_strand.items():
                        for v1_tm in v_dict_tm.values():
                            for v2_ivl in v1_tm:
                                info_tr_id.append(v2_ivl.transcript_model.trid)  # “info/ivls/tr_id“,
                                info_features_gene.append(
                                    v2_ivl.transcript_model.genename
                                )  # “info/ivls/features_gene“,
                                info_is_last3prime.append(v2_ivl.is_last_3prime)  # “info/ivls/is_last3prime“
                                info_is_intron.append(v2_ivl.kind == 105)  # “info/ivls/is_intron“,
                                info_start_end.append((v2_ivl.start, v2_ivl.end))  # “info/ivls/feture_start_end“
                                info_exino.append(v2_ivl.exin_no)  # “info/ivls/exino“
                                info_strandplus.append(
                                    v2_ivl.transcript_model.chromstrand[-1:] == "+"
                                )  # “info/ivls/strandplus“
                                info_chrm.append(v2_ivl.transcript_model.chromstrand[:-1])  # “info/ivls/chrm“

                    self.inv_tridstart2ix: dict[str, int] = {
                        f"{info_tr_id[i]}_{info_start_end[i][0]}": i for i in range(len(info_tr_id))
                    }
                    f.create_dataset(
                        "info/tr_id",
                        data=np.array(info_tr_id, dtype="S24"),
                        maxshape=(len(info_tr_id),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/features_gene",
                        data=np.array(info_features_gene, dtype="S15"),
                        maxshape=(len(info_features_gene),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_last3prime",
                        data=np.array(info_is_last3prime, dtype=bool),
                        maxshape=(len(info_is_last3prime),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_intron",
                        data=np.array(info_is_intron, dtype=bool),
                        maxshape=(len(info_is_intron),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/start_end",
                        data=np.array(info_start_end, dtype=np.int64),
                        maxshape=(len(info_start_end), 2),
                        chunks=(500, 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/exino",
                        data=np.array(info_exino, dtype=np.uint8),
                        maxshape=(len(info_exino),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/strandplus",
                        data=np.array(info_strandplus, dtype=np.bool),
                        maxshape=(len(info_strandplus),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/chrm",
                        data=np.array(info_chrm, dtype="S6"),
                        maxshape=(len(info_chrm),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )

                # cell_name = next(iter(molitems.keys())).split("$")[0]
                pos: DefaultDict[str, list[tuple[int, int]]] = defaultdict(list)
                mol: DefaultDict[str, list[int]] = defaultdict(list)
                ixs: DefaultDict[str, list[int]] = defaultdict(list)
                count_i: int = 0
                for mol_bc, molitem in molitems.items():
                    cell_name = mol_bc.split("$")[0]
                    with contextlib.suppress(StopIteration):
                        for match in next(iter(molitem.mappings_record.items()))[1]:
                            mol[cell_name].append(count_i)
                            pos[cell_name].append(match.segment)
                            ixs[cell_name].append(
                                self.inv_tridstart2ix[f"{match.feature.transcript_model.trid}_{match.feature.start}"]
                            )
                        count_i += 1
                # Do the last cell and close the file
                for cell_name in mol.keys():
                    posA = np.array(pos[cell_name], dtype=np.int32)
                    ixsA = np.array(ixs[cell_name], dtype=np.intp)
                    molA = np.array(mol[cell_name], dtype=np.uint32)
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/pos",
                        data=posA,
                        maxshape=posA.shape,
                        chunks=(min(500, posA.shape[0]), 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/ixs",
                        data=ixsA,
                        maxshape=ixsA.shape,
                        chunks=(min(500, ixsA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/mol",
                        data=molA,
                        maxshape=molA.shape,
                        chunks=(min(500, molA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                f.close()

        self.report_state += 1
        idx2bc = {v: k for k, v in bc2idx.items()}

        return dict_layers_columns, [idx2bc[i] for i in range(len(idx2bc))]

    def _count_cell_batch_stranded_discordant(
        self,
    ) -> tuple[dict[str, np.ndarray], list[str]]:
        """It performs molecule counting for the current batch of cells in the case of stranded method with discordant masking

        Returns
        -------
        dict_layers_columns: dict[str, np.ndarray]
            name_layer->np.ndarray of the batch
        idx2bc: list[str]
            list of barcodes

        NOTE This duplications of method is bad for code mantainance
        """
        molitems: DefaultDict[str, Molitem] = defaultdict(Molitem)
        # Sort similarly to what the sort linux command would do. (implemented using Read.__lt__)
        self.reads_to_count.sort()
        # NOTE: I could start by sorting the reads by chromosome, strand, position but for now let's see if it is fast without doing do

        repeats_reads_count = 0
        for r in self.reads_to_count:
            # Consider the correct strand
            ii = self.feature_indexes[f"{r.chrom}{r.strand}"]
            iir = self.feature_indexes[f"{r.chrom}{reverse(r.strand)}"]
            iim = self.mask_indexes[f"{r.chrom}{r.strand}"]
            iimr = self.mask_indexes[f"{r.chrom}{reverse(r.strand)}"]

            # Check if read is fully inside a masked region, in that case check if it is allowed in the reverse
            if iim.has_ivls_enclosing(r):
                repeats_reads_count += 1  # VERBOSE
                if not iimr.has_ivls_enclosing(r):
                    mappings_record = iir.find_overlapping_ivls(r)
                else:
                    continue
            else:
                # Look for overlap between the intervals and the read
                mappings_record = ii.find_overlapping_ivls(r)

            if len(mappings_record):
                # logger.debug("IN: Non empty mapping record")
                bcumi = f"{r.bc}${r.umi}"
                molitems[bcumi].add_mappings_record(mappings_record)
                # if len(molitems[bcumi].mappings_record):
                #     logger.debug("OUT: Non empty")
                # else:
                #     logger.debug("OUT: Empty")
        logger.debug(
            f"{repeats_reads_count} reads not considered because fully enclosed in repeat masked regions"
        )  # VERBOSE
        # initialize np.ndarray
        shape = (len(self.geneid2ix), len(self.cell_batch))

        dict_layers_columns: dict[str, np.ndarray] = {
            layer_name: sp.sparse.lil_array(shape, dtype=self.loom_numeric_dtype) for layer_name in self.logic.layers
        }
        bc2idx: dict[str, int] = dict(zip(self.cell_batch, range(len(self.cell_batch))))
        # After the whole file has been read, do the actual counting
        for bcumi, molitem in molitems.items():
            bc = bcumi.split("$")[0]  # extract the bc part from the bc+umi
            bcidx = bc2idx[bc]
            self.logic.count(molitem, bcidx, dict_layers_columns, self.geneid2ix)
            # NOTE I need to generalize this to any set of layers
            # before it was molitem.count(bcidx, spliced, unspliced, ambiguous, self.geneid2ix)

        if self.every_n_report and ((self.report_state % self.every_n_report) == 0):
            if self.kind_of_report == "p":
                import pickle

                first_cell_batch = next(iter(molitems.keys())).split("$")[0]
                if not (x := Path("pickle_dump")).exists():
                    x.mkdir()
                pickle.dump(
                    molitems,
                    open(f"pickle_dump/molitems_dump_{first_cell_batch}.pickle", "wb"),
                )
                pickle.dump(
                    self.reads_to_count,
                    open(f"pickle_dump/reads_to_count{first_cell_batch}.pickle", "wb"),
                )
            else:
                if not (x := Path(self.outputfolder).joinpath("dump")).exists():
                    x.mkdir()
                f = h5py.File(
                    Path(self.outputfolder).joinpath("dump", f"{self.sampleid}.hdf5")
                )  # From the docs: Read/write if exists, create otherwise (default)

                if "info/tr_id" not in f:
                    logger.warning(
                        "The hdf5 report is less accurate in reporting exactly all the information than the pickle."
                    )
                    info_tr_id = []
                    info_features_gene = []
                    info_is_last3prime = []
                    info_is_intron = []
                    info_start_end = []
                    info_exino = []
                    info_strandplus = []
                    info_chrm = []
                    for _, v_dict_tm in self.annotations_by_chrm_strand.items():
                        for v1_tm in v_dict_tm.values():
                            for v2_ivl in v1_tm:
                                info_tr_id.append(v2_ivl.transcript_model.trid)  # “info/ivls/tr_id“,
                                info_features_gene.append(
                                    v2_ivl.transcript_model.genename
                                )  # “info/ivls/features_gene“,
                                info_is_last3prime.append(v2_ivl.is_last_3prime)  # “info/ivls/is_last3prime“
                                info_is_intron.append(v2_ivl.kind == 105)  # “info/ivls/is_intron“,
                                info_start_end.append((v2_ivl.start, v2_ivl.end))  # “info/ivls/feture_start_end“
                                info_exino.append(v2_ivl.exin_no)  # “info/ivls/exino“
                                info_strandplus.append(
                                    v2_ivl.transcript_model.chromstrand[-1:] == "+"
                                )  # “info/ivls/strandplus“
                                info_chrm.append(v2_ivl.transcript_model.chromstrand[:-1])  # “info/ivls/chrm“

                    self.inv_tridstart2ix: dict[str, int] = {
                        f"{info_tr_id[i]}_{info_start_end[i][0]}": i for i in range(len(info_tr_id))
                    }
                    f.create_dataset(
                        "info/tr_id",
                        data=np.array(info_tr_id, dtype="S24"),
                        maxshape=(len(info_tr_id),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/features_gene",
                        data=np.array(info_features_gene, dtype="S15"),
                        maxshape=(len(info_features_gene),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_last3prime",
                        data=np.array(info_is_last3prime, dtype=bool),
                        maxshape=(len(info_is_last3prime),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_intron",
                        data=np.array(info_is_intron, dtype=bool),
                        maxshape=(len(info_is_intron),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/start_end",
                        data=np.array(info_start_end, dtype=np.int64),
                        maxshape=(len(info_start_end), 2),
                        chunks=(500, 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/exino",
                        data=np.array(info_exino, dtype=np.uint8),
                        maxshape=(len(info_exino),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/strandplus",
                        data=np.array(info_strandplus, dtype=np.bool),
                        maxshape=(len(info_strandplus),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/chrm",
                        data=np.array(info_chrm, dtype="S6"),
                        maxshape=(len(info_chrm),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )

                # cell_name = next(iter(molitems.keys())).split("$")[0]
                pos: DefaultDict[str, list[tuple[int, int]]] = defaultdict(list)
                mol: DefaultDict[str, list[int]] = defaultdict(list)
                ixs: DefaultDict[str, list[int]] = defaultdict(list)
                count_i: int = 0
                for mol_bc, molitem in molitems.items():
                    cell_name = mol_bc.split("$")[0]
                    with contextlib.suppress(StopIteration):
                        for match in next(iter(molitem.mappings_record.items()))[1]:
                            mol[cell_name].append(count_i)
                            pos[cell_name].append(match.segment)
                            ixs[cell_name].append(
                                self.inv_tridstart2ix[f"{match.feature.transcript_model.trid}_{match.feature.start}"]
                            )
                        count_i += 1
                # Do the last cell and close the file
                for cell_name in mol.keys():
                    posA = np.array(pos[cell_name], dtype=np.int32)
                    ixsA = np.array(ixs[cell_name], dtype=np.intp)
                    molA = np.array(mol[cell_name], dtype=np.uint32)
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/pos",
                        data=posA,
                        maxshape=posA.shape,
                        chunks=(min(500, posA.shape[0]), 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/ixs",
                        data=ixsA,
                        maxshape=ixsA.shape,
                        chunks=(min(500, ixsA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/mol",
                        data=molA,
                        maxshape=molA.shape,
                        chunks=(min(500, molA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                f.close()

        self.report_state += 1
        idx2bc = {v: k for k, v in bc2idx.items()}

        return dict_layers_columns, [idx2bc[i] for i in range(len(idx2bc))]

    def _count_cell_batch_non_stranded(self) -> tuple[dict[str, np.ndarray], list[str]]:
        """It performs molecule counting for the current batch of cells in the case of non stranded method

        Returns
        -------
        dict_layers_columns: dict[str, np.ndarray]
            name_layer->np.ndarray of the batch
        idx2bc: list[str]
            list of barcodes
        """
        molitems: DefaultDict[str, Molitem] = defaultdict(Molitem)
        # Sort similarly to what the sort linux command would do. (implemented using Read.__lt__)
        # NOTE NOTE!!!! Here I changed the way to sort because it was using the strand causing to skip a lot of reads in SmartSeq2
        self.reads_to_count.sort()
        # NOTE: I could start by sorting the reads by chromosome, strand, position but for now let's see if it is fast without doing do

        repeats_reads_count = 0
        plus_reads_count = 0
        minus_reads_count = 0
        both_reads_count = 0
        for r in self.reads_to_count:
            # Consider the correct strand
            ii = self.feature_indexes[f"{r.chrom}{r.strand}"]
            iir = self.feature_indexes[f"{r.chrom}{reverse(r.strand)}"]
            iim = self.mask_indexes[f"{r.chrom}{r.strand}"]
            iimr = self.mask_indexes[f"{r.chrom}{reverse(r.strand)}"]

            # Check if read is fully inside a masked region, in that case skip it
            if iim.has_ivls_enclosing(r) or iimr.has_ivls_enclosing(r):
                repeats_reads_count += 1  # VERBOSE
                continue

            # Look for overlap between the intervals and the read
            mappings_record = ii.find_overlapping_ivls(r)
            if len(mappings_record):
                bcumi = f"{r.bc}${r.umi}"
                molitems[bcumi].add_mappings_record(mappings_record)
                if r.strand == "+":
                    plus_reads_count += 1
                else:
                    minus_reads_count += 1

            mappings_record_r = iir.find_overlapping_ivls(r)
            if len(mappings_record_r):
                bcumi = f"{r.bc}${r.umi}"
                molitems[bcumi].add_mappings_record(mappings_record_r)
                if r.strand == "-":
                    plus_reads_count += 1
                else:
                    minus_reads_count += 1

            if len(mappings_record) and len(mappings_record_r):
                both_reads_count += 1

        logger.debug(f"{repeats_reads_count} reads in repeat masked regions")  # VERBOSE
        logger.debug(f"{plus_reads_count} reads overlapping with features on plus strand")  # VERBOSE
        logger.debug(f"{minus_reads_count} reads overlapping with features on minus strand")  # VERBOSE
        logger.debug(f"{both_reads_count} reads overlapping with features on both strands")  # VERBOSE
        # initialize np.ndarray
        shape = (len(self.geneid2ix), len(self.cell_batch))

        dict_layers_columns: dict[str, np.ndarray] = {
            layer_name: sp.sparse.lil_array(shape, dtype=self.loom_numeric_dtype) for layer_name in self.logic.layers
        }
        bc2idx: dict[str, int] = dict(zip(self.cell_batch, range(len(self.cell_batch))))
        # After the whole file has been read, do the actual counting
        for bcumi, molitem in molitems.items():
            bc = bcumi.split("$")[0]  # extract the bc part from the bc+umi
            bcidx = bc2idx[bc]
            self.logic.count(molitem, bcidx, dict_layers_columns, self.geneid2ix)
            # NOTE I need to generalize this to any set of layers
            # before it was molitem.count(bcidx, spliced, unspliced, ambiguous, self.geneid2ix)

        if self.every_n_report and ((self.report_state % self.every_n_report) == 0):
            if self.kind_of_report == "p":
                import pickle

                first_cell_batch = next(iter(molitems.keys())).split("$")[0]
                if not (x := Path("pickle_dump")).exists():
                    x.mkdir()
                pickle.dump(
                    molitems,
                    open(f"pickle_dump/molitems_dump_{first_cell_batch}.pickle", "wb"),
                )
                pickle.dump(
                    self.reads_to_count,
                    open(f"pickle_dump/reads_to_count{first_cell_batch}.pickle", "wb"),
                )
            else:
                if not (x := Path(self.outputfolder).joinpath("dump")).exists():
                    x.mkdir()
                f = h5py.File(
                    Path(self.outputfolder).joinpath("dump", f"{self.sampleid}.hdf5")
                )  # From the docs: Read/write if exists, create otherwise (default)

                if "info/tr_id" not in f:
                    logger.warning(
                        "The hdf5 report is less accurate than the pickle in the completeness of the info it is reporting."
                    )
                    info_tr_id = []
                    info_features_gene = []
                    info_is_last3prime = []
                    info_is_intron = []
                    info_start_end = []
                    info_exino = []
                    info_strandplus = []
                    info_chrm = []
                    for _k, v_dict_tm in self.annotations_by_chrm_strand.items():
                        for v1_tm in v_dict_tm.values():
                            for v2_ivl in v1_tm:
                                info_tr_id.append(v2_ivl.transcript_model.trid)  # “info/ivls/tr_id“,
                                info_features_gene.append(
                                    v2_ivl.transcript_model.genename
                                )  # “info/ivls/features_gene“,
                                info_is_last3prime.append(v2_ivl.is_last_3prime)  # “info/ivls/is_last3prime“
                                info_is_intron.append(v2_ivl.kind == 105)  # “info/ivls/is_intron“,
                                info_start_end.append((v2_ivl.start, v2_ivl.end))  # “info/ivls/feture_start_end“
                                info_exino.append(v2_ivl.exin_no)  # “info/ivls/exino“
                                info_strandplus.append(
                                    v2_ivl.transcript_model.chromstrand[-1:] == "+"
                                )  # “info/ivls/strandplus“
                                info_chrm.append(v2_ivl.transcript_model.chromstrand[:-1])  # “info/ivls/chrm“

                    self.inv_tridstart2ix: dict[str, int] = {
                        f"{info_tr_id[i]}_{info_start_end[i][0]}": i for i in range(len(info_tr_id))
                    }
                    f.create_dataset(
                        "info/tr_id",
                        data=np.array(info_tr_id, dtype="S24"),
                        maxshape=(len(info_tr_id),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/features_gene",
                        data=np.array(info_features_gene, dtype="S15"),
                        maxshape=(len(info_features_gene),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_last3prime",
                        data=np.array(info_is_last3prime, dtype=bool),
                        maxshape=(len(info_is_last3prime),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/is_intron",
                        data=np.array(info_is_intron, dtype=bool),
                        maxshape=(len(info_is_intron),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/start_end",
                        data=np.array(info_start_end, dtype=np.int64),
                        maxshape=(len(info_start_end), 2),
                        chunks=(500, 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/exino",
                        data=np.array(info_exino, dtype=np.uint8),
                        maxshape=(len(info_exino),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/strandplus",
                        data=np.array(info_strandplus, dtype=np.bool),
                        maxshape=(len(info_strandplus),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        "info/chrm",
                        data=np.array(info_chrm, dtype="S6"),
                        maxshape=(len(info_chrm),),
                        chunks=(500,),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )

                # cell_name = next(iter(molitems.keys())).split("$")[0]
                pos: DefaultDict[str, list[tuple[int, int]]] = defaultdict(list)
                mol: DefaultDict[str, list[int]] = defaultdict(list)
                ixs: DefaultDict[str, list[int]] = defaultdict(list)
                count_i: int = 0
                for mol_bc, molitem in molitems.items():
                    cell_name = mol_bc.split("$")[0]
                    with contextlib.suppress(StopIteration):
                        for match in next(iter(molitem.mappings_record.items()))[1]:
                            mol[cell_name].append(count_i)
                            pos[cell_name].append(match.segment)
                            ixs[cell_name].append(
                                self.inv_tridstart2ix[f"{match.feature.transcript_model.trid}_{match.feature.start}"]
                            )
                        count_i += 1
                # Do the last cell and close the file
                for cell_name in mol.keys():
                    posA = np.array(pos[cell_name], dtype=np.int32)
                    ixsA = np.array(ixs[cell_name], dtype=np.intp)
                    molA = np.array(mol[cell_name], dtype=np.uint32)
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/pos",
                        data=posA,
                        maxshape=posA.shape,
                        chunks=(min(500, posA.shape[0]), 2),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/ixs",
                        data=ixsA,
                        maxshape=ixsA.shape,
                        chunks=(min(500, ixsA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                    f.create_dataset(
                        f"cells/{self.sampleid}_{cell_name}/mol",
                        data=molA,
                        maxshape=molA.shape,
                        chunks=(min(500, molA.shape[0]),),
                        compression="gzip",
                        shuffle=False,
                        compression_opts=4,
                    )
                f.close()

        self.report_state += 1
        idx2bc = {v: k for k, v in bc2idx.items()}

        return dict_layers_columns, [idx2bc[i] for i in range(len(idx2bc))]

    def pcount(
        self,
        samfile: str,
        cell_batch_size: int = 50,
        molecules_report: bool = False,
        n_processes: int = 4,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        """Do the counting of molecules in parallel using multiprocessing"""

        raise NotImplementedError("Implement this using multiprocessiong")

    def pcount_cell_batch(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
        """It performs molecule counting for the current batch of cells"""
        raise NotImplementedError("This will be a used by .pcount")


def reverse(strand: str) -> str:
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"
    else:
        raise ValueError(f"Unknown strand {strand}")
