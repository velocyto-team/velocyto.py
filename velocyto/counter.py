import re
import sys
import gzip
import logging
from typing import *
from collections import defaultdict
from itertools import chain
from collections import OrderedDict
import velocyto as vcy
import pysam
import numpy as np
import os
import sys


class ExInCounter:
    """ Main class to do the counting of introns and exons """
    def __init__(self, logic: vcy.Logic, valid_bcset: Set[str]=None) -> None:
        """ valid_bcs2idx is a dict mapping valid cell barcodes to unique arbitrary numeric indexes used in the output arrays """
        self.logic = logic()
        # NOTE: maybe there shoulb be a self.logic.validate() step at the end of init
        if valid_bcset is None:
            self.valid_bcset: Set[str] = set()
            self.filter_mode = False
            self.counter = 0
        else:
            self.valid_bcset = valid_bcset
            self.filter_mode = True
        self.annotations_by_chrm_strand: Dict[str, Dict[str, vcy.TranscriptModel]] = {}
        self.mask_ivls_by_chromstrand = defaultdict(list)  # type: Dict[str, List]
        self.geneid2ix: Dict[str, int] = {}
        self.genes: Dict[str, vcy.GeneInfo] = {}
        # NOTE: by using a default dict and not logging access to keys that do not exist, we might miss bugs!!!
        self.test_flag = None

    @property
    def bclen(self) -> int:
        try:
            return len(next(iter(self.valid_bcset)))
        except StopIteration:
            return None
    
    @staticmethod
    def parse_cigar_tuple(cigartuples: List[Tuple], pos: int) -> Tuple[List[Tuple[int, int]], bool, int, int]:
        segments = []
        hole_to_remove = []
        ref_skip = False
        clip5 = clip3 = 0
        p = pos
        for operation_id, length in cigartuples:
            if operation_id == 0:  # vcy.CIGAR[operation_id] == "BAM_CMATCH"
                segments.append((p, p + length - 1))
                p += length
            elif operation_id == 3:  # A splice || vcy.CIGAR[operation_id] == 'BAM_CREF_SKIP'
                ref_skip = True
                p += length
            elif operation_id == 2:  # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
                if length <= vcy.PATCH_INDELS:
                    hole_to_remove.append(len(segments) - 1)
                p += length
            elif operation_id == 4:  # bases at 5' or 3' are NOT part of the alignment || vcy.CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
                if p == pos:
                    clip5 = length  # At start of alignment
                else:
                    clip3 = length  # Must be at end of alignment vcy.CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
                p += length
            elif operation_id == 1:  # An insertion BAM_CINS
                if length <= vcy.PATCH_INDELS:
                    hole_to_remove.append(len(segments) - 1)
                # else do nothing
                # NOTE: maybe we should make so that the reads get discarded
            elif operation_id == 5:  # BAM_CHARD_CLIP
                logging.warn("Hard clip was encountered! All mapping are assumed soft clipped")

        # Merge segments separated by small insertions and deletions
        for a, b in enumerate(hole_to_remove):
            segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])

        return segments, ref_skip, clip5, clip3

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

    def iter_alignments(self, samfile: str, unique: bool=True, yield_line: bool=False) -> Iterable:
        """Iterates over the bam/sam file and yield Read objects

        Arguments
        ---------
        samfile: str
            path to the sam file
        unique: bool
            yield only unique allignments
        yield_line: bool
            whether to yield the raw sam line

        Returns
        -------
        yields vcy.Read for each valid line of the bam file
        or a Tuple (vcy.Read, sam_line) if ``yield_line==True``
        """
        # No peeking here, this will happen a layer above, and only once  on the not sorted file! before it was self.peek(samfile, lines=10)
        logging.debug(f"Reading {samfile}")
        fin = pysam.AlignmentFile(samfile)  # type: pysam.AlignmentFile
        for i, read in enumerate(fin):
            if i % 10000000 == 0:
                logging.debug(f"Read first {i // 1000000} million reads")
            if read.is_unmapped:
                continue
            # If unique parameter is set to True, skip not unique allignments
            if unique and read.get_tag("NH") != 1:
                continue
            try:
                bc = read.get_tag(self.cellbarcode_str).split("-")[0]  # NOTE: this rstrip is relevant only for cellranger, should not cause trouble in Dropseq
                umi = read.get_tag(self.umibarcode_str)
            except KeyError:
                continue  # NOTE: Here errors could go unnoticed
            if bc not in self.valid_bcset:
                if self.filter_mode:
                    continue
                else:
                    self.valid_bcset.add(bc)
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
            pos = read.reference_start + 1  # reads in pysam are always 0-based, but 1-based is more convenient to wor with in bioinformatics
            segments, ref_skipped, clip5, clip3 = self.parse_cigar_tuple(read.cigartuples, pos)
            if segments == []:
                logging.debug("No segments in read:%s" % read.qname)

            read_object = vcy.Read(bc, umi, chrom, strand, pos, segments, clip5, clip3, ref_skipped)
            if yield_line:
                if read_object.span > 3000000:  # Longest locus existing
                    logging.warn(f"Trashing read, too long span\n{read.tostring(fin)}")
                else:
                    yield read_object, read.tostring(fin)
            else:
                if read_object.span > 3000000:  # Longest locus existing
                    logging.warn(f"Trashing read, too long span\n{read.tostring(fin)}")
                else:
                    yield read_object
        fin.close()

    def read_repeats(self, gtf_file: str, tolerance: int=5) -> Dict[str, List[vcy.Feature]]:
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
        mask_ivls_by_chromstrand: Dict[str, List[vcy.Feature]]
            A dictionary key: chromosome+strand value: list of features (repeat intervals)
            (The reference is returned but an internal attribure self.self.masked_by_chrm_strand is kept)
        
        """
        # Example code to sort the gtf file
        # sorted_filename = gtffile.split(".")[0] + "_sorted.gtf"
        # logging.debug(f"Sorting by `sort -k1,1 -k7,7 -k4,4n {gtffile} > {sorted_filename}`")
        # with open(sorted_filename, "w") as f:
        #     p1 = subprocess.run(["sort", "-k1,1", "-k7,7", "-k4,4n", gtffile],
        #                         stdout=f)
            
        # Parse arguments
        logging.debug(f'Reading {gtf_file}, the file is assumed sorted (if is not interupt the process and sort it by running: ``sort -k1,1 -k7,7 -k4,4n``)')
        
        # Read up skipping headers up to the first valid entry
        repeat_ivls_list: List[vcy.Feature] = []

        # fin = open(gtf_file)
        gtf_lines = [line for line in open(gtf_file) if not line.startswith('#')]

        def sorting_key(entry: str) -> Tuple[str, bool, int, str]:
            """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
            x = entry.split("\t")
            return (x[0], x[6] == "+", int(x[3]), entry)  # The last element of the touple corresponds to the `last resort comparison`
        
        gtf_lines = sorted(gtf_lines, key=sorting_key)
        ###

        line = gtf_lines.pop(0)
        fields = line.rstrip().split('\t')
        chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
        # Removing chr from the chromosome name to uniform different formats of gtf files, taht might or might not have the prefix "chr"
        if chrom[:3].lower() == "chr":
            chrom = chrom[3:]
        start = int(start_str)
        end = int(end_str)
        chromstrand: str = chrom + strand

        # Set this tu the current entry
        curr_chrom: str = chrom
        curr_feature_class: str = feature_class
        curr_feature_type: str = feature_type
        curr_start: int = start
        curr_end: int = end
        curr_n: int = 1
        curr_strand: str = strand
        curr_tags: str = tags
        curr_chromstrand: str = chromstrand

        for line in gtf_lines:
            
            fields = line.rstrip().split('\t')
            chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
            # Removing chr from the chromosome name to uniform different formats of gtf files, taht might or might not have the prefix "chr"
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
                repeat_ivls_list.append(vcy.Feature(start=curr_start, end=curr_end, kind=ord("r"), exin_no=curr_n))
                # NOTE exin_no is being used to store curr_n, and kind could be used to store tag for some extra memory cost
                
                # Remember the newly parsed interval
                curr_start = start
                curr_end = end
                curr_n = 1  # number of original repeat regions included

                # this is extra information that right now is not stored in vcy.Feature (not to waste memory) since is not used
                curr_feature_class = feature_class
                curr_feature_type = feature_type
                curr_tags = tags
            else:
                # Extend the masked interval!
                curr_end = end
                curr_n += 1
                # this is extra information that right now is not saved (not to waste memory) since is not used
                gap = start - curr_end
                curr_tags = f"{curr_tags} gap {gap}; {tags}" if gap > 0 else curr_tags + tags

        n = 0
        for chromstrand, feature_list in self.mask_ivls_by_chromstrand.items():
            feature_list.sort()  # relies on the __lt__ method of vcy.Feature
            n += len(feature_list)

        logging.debug(f'Processed masked annotation .gtf and generated {n} intervals to mask!')

        return self.mask_ivls_by_chromstrand

    def assign_indexes_to_genes(self, features: Dict[str, vcy.TranscriptModel]) -> None:
        """Assign to each newly encoutered genes an unique indexes corresponding to the output matrix column ix
        """
        logging.debug("Assigning indexes to genes")
        for name, trmodel in features.items():
            if trmodel.geneid in self.geneid2ix:
                if self.genes[trmodel.geneid].start > trmodel.start:
                    self.genes[trmodel.geneid].start = trmodel.start
                if self.genes[trmodel.geneid].end < trmodel.end:
                    self.genes[trmodel.geneid].end = trmodel.end
            else:
                self.geneid2ix[trmodel.geneid] = len(self.geneid2ix)
                self.genes[trmodel.geneid] = vcy.GeneInfo(trmodel.genename, trmodel.geneid, trmodel.chromstrand, trmodel.start, trmodel.end)

    def read_transcriptmodels(self, gtf_file: str) -> Dict[str, Dict[str, vcy.TranscriptModel]]:
        """Reads transcript models from a sorted .gtf file

        Arguments
        ---------
        gtf_file: str
            Path to the sorted gtf file

        Returns
        -------
        annotations_by_chrm_strand: Dict[str, List[vcy.TrancriptModel]]
            A dictionary key: chromosome+strand value: list of trascript models
            (The reference is returned but an internal attribure self.annotations_by_chrm_strand is kept)

        There will exist an object vcy.Features for the same exon appearing in a different vcy.TranscriptModel. (his is desired)
        """

        # Define the regexes for the parsing
        regex_trid = re.compile('transcript_id "([^"]+)"')
        regex_trname = re.compile('transcript_name "([^"]+)"')
        regex_geneid = re.compile('gene_id "([^"]+)"')
        regex_genename = re.compile('gene_name "([^"]+)"')
        regex_exonno = re.compile('exon_number "([^"]+)"')

        # Initialize containers
        # headerlines: List[str] = []
        curr_chromstrand = None
        features: Dict[str, vcy.TranscriptModel] = OrderedDict()

        gtf_lines = [line for line in open(gtf_file) if not line.startswith('#')]

        def sorting_key(entry: str) -> Tuple[str, bool, int, str]:
            """This sorting strategy is equivalent to sort -k1,1 -k7,7 -k4,4n"""
            x = entry.split("\t")
            return (x[0], x[6] == "+", int(x[3]), entry)  # The last element of the touple corresponds to the `last resort comparison`
        
        gtf_lines = sorted(gtf_lines, key=sorting_key)
        # Loop throug gtf file (assumes it is ordered)
        for nth_line, line in enumerate(gtf_lines):
            # Deal with headers
            # if line.startswith('#'):
            #     headerlines.append(line)
            #     continue
                
            fields = line.rstrip().split('\t')
            chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields

            # Deal with possible incongruences between .gtf and bam file chromosome naming
            # We keep the name of the chromosome removing the chr from the name and removing part after the . if present
            if "chr" in chrom[:4]:
                chrom = chrom[3:]  # NOTE before it was chrom[3:].split(".")[0]
            else:
                pass  # NOTE before it was chrom.split(".")[0]

            # A new chromosome/strand encountered
            if chrom + strand != curr_chromstrand:
                if curr_chromstrand is not None:  # Every time with exception with first and the last chromosome
                    if not chrom + strand not in self.annotations_by_chrm_strand:
                        # NOTE this is not enough as a check but it will detect with few checks if file is not sorted at all
                        raise IOError(f"Genome annotation gtf file is not sorted correctly! Run the following command:\nsort -k1,1 -k7,7 -k4,4n -o [GTF_OUTFILE] [GTF_INFILE]")
                    else:
                        logging.debug(f"Done with {curr_chromstrand} [line {nth_line-1}]")
                    self.assign_indexes_to_genes(features)
                    self.annotations_by_chrm_strand[curr_chromstrand] = features
                    logging.debug(f"Seen {len(self.geneid2ix)} genes until now")
                features = OrderedDict()
                logging.debug(f"Parsing Chromosome {chrom} strand {strand} [line {nth_line}]")
                curr_chromstrand = chrom + strand
                
            if feature_type in ("exon"):
                trid = regex_trid.search(tags).group(1)
                trname = regex_trname.search(tags).group(1)
                geneid = regex_geneid.search(tags).group(1)
                genename = regex_genename.search(tags).group(1)
                exonno = regex_exonno.search(tags).group(1)
                start = int(start_str)
                end = int(end_str)
                chromstrand = chrom + strand
                
                try:
                    features[trid].append_exon(vcy.Feature(start=start, end=end, kind=ord("e"), exin_no=exonno))
                except KeyError:
                    features[trid] = vcy.TranscriptModel(trid=trid, trname=trname, geneid=geneid, genename=genename, chromstrand=chromstrand)
                    features[trid].append_exon(vcy.Feature(start=start, end=end, kind=ord("e"), exin_no=exonno))

        # Do it for the last chromosome
        self.assign_indexes_to_genes(features)
        self.annotations_by_chrm_strand[curr_chromstrand] = features
        logging.debug(f"Done with {curr_chromstrand} [line {nth_line-1}]")
    
        logging.debug(f"Fixing corner cases of transcript models containg intron longer than {vcy.LONGEST_INTRON_ALLOWED//1000}Kbp")
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

    def mark_up_introns(self, bamfile: str, multimap: bool) -> None:
        """ Mark up introns that have reads across exon-intron junctions
        
        Arguments
        ---------
        bamfile: str
            path to the bam or sam file to markup
        logic: vcy.Logic
            The logic object to use, changes in different techniques / levels of strictness
            NOTE: Right now it is not used

        Returns
        -------
        Nothing it just add to validation to the vcy.Features

        Note
        ----
        Situations not considered:
        # If an the exon is so short that is possible to get both exonA-exonB junction and exonB-intronB boundary in the same read
        """

        # VERBOSE: dump_list = []
        # Read the file
        currchrom = ""
        set_chromosomes_seen: Set[str] = set()
        for r in self.iter_alignments(bamfile, unique=not multimap):
            # Don't consider spliced reads (exonic) in this step
            # NOTE Can the exon be so short that we get splicing and exon-intron boundary
            if r.is_spliced:
                continue

            # When the chromosome mapping of the read changes, change interval index.
            if r.chrom != currchrom:
                if r.chrom in set_chromosomes_seen:
                    raise IOError(f"Input .bam file should be chromosome-sorted. (Hint: use `samtools sort {bamfile}`)")
                set_chromosomes_seen.add(r.chrom)
                logging.debug(f"Marking up chromosome {r.chrom}")
                currchrom = r.chrom
                if currchrom + '+' not in self.annotations_by_chrm_strand:
                    logging.warn(f"The .bam file refers to a chromosome '{currchrom}+' not present in the annotation (.gtf) file")
                    sorted_features_f: List[vcy.Feature] = []
                else:
                    sorted_features_f = sorted(chain.from_iterable(self.annotations_by_chrm_strand[currchrom + '+'].values()))
                if currchrom + '-' not in self.annotations_by_chrm_strand:
                    logging.warn(f"The .bam file refers to a chromosome '{currchrom}-' not present in the annotation (.gtf) file")
                    sorted_features_r: List[vcy.Feature] = []
                else:
                    sorted_features_r = sorted(chain.from_iterable(self.annotations_by_chrm_strand[currchrom + '-'].values()))
                iif = vcy.FeatureIndex(sorted_features_f)
                iir = vcy.FeatureIndex(sorted_features_r)
            
            # Consider the correct strand
            ii = iif if r.strand == '+' else iir

            # VERBOSE: # Look for overlap between the intervals and the read
            # VERBOSE: dump_list += ii.mark_overlapping_ivls(r)
            ii.mark_overlapping_ivls(r)

        # VERBOSE: import pickle
        # VERBOSE: pickle.dump(dump_list, open("dump_mark_overlapping_ivls.pickle", "wb"))

    def count(self, samfile: str, multimap: bool, cell_batch_size: int=100, molecules_report: bool=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """ Do the counting of molecules
        
        Arguments
        ---------
        samfile: str
            path to the bam or sam file to markup
        cell_batch_size: int, default = 50
            it defines whether to require or not exon-intron spanning read to consider an intron valid.
        molecules_report: bool, default=False
            whether to output a file that reports the intervals hit for each molecule counted
        
        Returns
        -------
        list_spliced_arrays, list_unspliced_arrays, list_ambiguous_arrays, cell_bcs_order

        Note
        ----
        The memory footprint could be reduced allocating scipy sparse arrays
        
        """

        # Initialize variables and containers
        self.cell_batch: Set[str] = set()
        self.reads_to_count: List[vcy.Read] = []
        # self.cells_since_last_count = 0
        # Analysis is cell wise so the Feature Index swapping is happening often and it is worth to preload everything in memory
        self.feature_indexes: DefaultDict[str, vcy.FeatureIndex] = defaultdict(vcy.FeatureIndex)
        for chromstrand_key, annotions_ordered_dict in self.annotations_by_chrm_strand.items():
            self.feature_indexes[chromstrand_key] = vcy.FeatureIndex(sorted(chain.from_iterable(annotions_ordered_dict.values())))

        self.mask_indexes: DefaultDict[str, vcy.FeatureIndex] = defaultdict(vcy.FeatureIndex)
        for chromstrand_key, annotions_list in self.mask_ivls_by_chromstrand.items():
            self.mask_indexes[chromstrand_key] = vcy.FeatureIndex(annotions_list)  # This suould be sorted

        # Before counting, report how many features where validated
        logging.debug(f"Summarizing the results of intron validation.")
        n_is_intron = 0
        n_is_intron_valid = 0
        unique_valid = set()
        for chromstrand_key, feature_index in self.feature_indexes.items():
            for ivl in feature_index.ivls:
                if ivl.kind == ord("i"):
                    n_is_intron += 1
                if ivl.is_validated:
                    n_is_intron_valid += 1
                    unique_valid.add((ivl.start, ivl.end))
        logging.debug(f"Validated {n_is_intron_valid} introns (of which unique intervals {len(unique_valid)}) out of {n_is_intron} total possible introns (considering each possible transcript models).")

        cell_bcs_order: List[str] = []
        list_spliced_arrays: List[np.ndarray] = []
        list_unspliced_arrays: List[np.ndarray] = []
        list_ambiguous_arrays: List[np.ndarray] = []
        nth = 0
        # Loop through the aligment of the samfile
        for r in self.iter_alignments(samfile, unique=not multimap):
            if len(self.cell_batch) == cell_batch_size and r.bc not in self.cell_batch:
                # Perfrom the molecule counting
                nth += 1
                logging.debug(f"Counting for batch {nth}, containing {len(self.cell_batch)} cells and {len(self.reads_to_count)} reads")
                spliced, unspliced, ambiguous, list_bcs = self.count_cell_batch(molecules_report and (nth % 10 == 0))
                cell_bcs_order += list_bcs
                list_spliced_arrays.append(spliced)
                list_unspliced_arrays.append(unspliced)
                list_ambiguous_arrays.append(ambiguous)
                self.cell_batch = set()
                # Drop the counted reads (If there are no other reference to it) and reset the indexes to 0
                self.reads_to_count = []
                for chromstrand_key, annotions_ordered_dict in self.annotations_by_chrm_strand.items():
                    self.feature_indexes[chromstrand_key].reset()
                for chromstrand_key, annotions_list in self.mask_ivls_by_chromstrand.items():
                    self.mask_indexes[chromstrand_key].reset()

            self.cell_batch.add(r.bc)
            self.reads_to_count.append(r)
        logging.debug(f"Counting molecule for last batch of {len(self.cell_batch)}")
        spliced, unspliced, ambiguous, list_bcs = self.count_cell_batch()
        cell_bcs_order += list_bcs
        list_spliced_arrays.append(spliced)
        list_unspliced_arrays.append(unspliced)
        list_ambiguous_arrays.append(ambiguous)
        self.cell_batch = set()
        self.reads_to_count = []
        logging.debug(f"Counting done!")
        return list_spliced_arrays, list_unspliced_arrays, list_ambiguous_arrays, cell_bcs_order

    def count_cell_batch(self, molecules_report: bool=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """It performs molecule counting for the current batch of cells

        Returns
        -------
        spliced: np.ndarray
        unspliced: np.ndarray
        ambiguous: np.ndarray
        idx2bc: List[str]
            list of barcodes
        """
        molitems: DefaultDict[str, vcy.Molitem] = defaultdict(vcy.Molitem)
        # Sort similarly to what the sort linux command would do. (implemented using Read.__lt__)
        self.reads_to_count.sort()
        # NOTE: I could start by sorting the reads by chromosome, strand, position but for now let's see if it is fast without doing do

        repeats_reads_count = 0
        for r in self.reads_to_count:
            # Consider the correct strand
            ii = self.feature_indexes[r.chrom + r.strand]
            iim = self.mask_indexes[r.chrom + r.strand]

            # Check if read is fully inside a repeat, in that case skip it
            if iim.has_ivls_enclosing(r):
                repeats_reads_count += 1  # VERBOSE
                continue

            # Look for overlap between the intervals and the read
            mappings_record = ii.find_overlapping_ivls(r)
            if len(mappings_record):
                # logging.debug("IN: Non empty mapping record")
                bcumi = r.bc + r.umi
                molitems[bcumi].add_mappings_record(mappings_record)
                # if len(molitems[bcumi].mappings_record):
                #     logging.debug("OUT: Non empty")
                # else:
                #     logging.debug("OUT: Empty")
        logging.debug(f"{repeats_reads_count} reads not considered because fully enclosed in repeat masked regions")  # VERBOSE
        # initialize np.ndarray
        shape = (len(self.geneid2ix), len(self.cell_batch))
        spliced = np.zeros(shape, dtype=vcy.LOOM_NUMERIC_DTYPE)
        unspliced = np.zeros(shape, dtype=vcy.LOOM_NUMERIC_DTYPE)
        ambiguous = np.zeros(shape, dtype=vcy.LOOM_NUMERIC_DTYPE)
        bc2idx: Dict[str, int] = dict(zip(self.cell_batch, range(len(self.cell_batch))))

        # After the whole file has been read, do the actual counting
        for bcumi, molitem in molitems.items():
            bc = bcumi[:self.bclen]  # extract the bc part from the bc+umi
            bcidx = bc2idx[bc]
            self.logic.count(molitem, bcidx, spliced, unspliced, ambiguous, self.geneid2ix)
            # NOTE I need to generalize this to any set of layers
            # before it was molitem.count(bcidx, spliced, unspliced, ambiguous, self.geneid2ix)

        if molecules_report:
            import pickle
            first_cell_batch = next(iter(molitems.keys()))[:self.bclen]
            if not os.path.exists("pickle_dump"):
                os.makedirs("pickle_dump")
            pickle.dump(molitems, open(f"pickle_dump/molitems_dump_{first_cell_batch}.pickle", "wb"))
            pickle.dump(self.reads_to_count, open(f"pickle_dump/reads_to_count{first_cell_batch}.pickle", "wb"))

        idx2bc = {v: k for k, v in bc2idx.items()}
        
        return spliced, unspliced, ambiguous, [idx2bc[i] for i in range(len(idx2bc))]

    def pcount(self, samfile: str, cell_batch_size: int=50, molecules_report: bool=False, n_processes: int=4) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """ Do the counting of molecules in parallel using multiprocessing
        """

        raise NotImplementedError("Implement this using multiprocessiong")

    def pcount_cell_batch(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
        """It performs molecule counting for the current batch of cells
        """
        raise NotImplementedError("This will be a used by .pcount")
