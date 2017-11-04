import logging
import copy
import os
import click
import re
import sys
import gzip
import bisect
from collections import defaultdict
from typing import *
import velocyto as vcy

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


@click.command(short_help="Transform a genome annotation .gtf file into a intervals .txt file.")
@click.option('--dogenes/--no-dogenes',
              default=True,
              help="whether to process gene models",
              show_default=True)
@click.option('--dotrs/--no-dotrs',
              default=True,
              help="whether to process transcript models",
              show_default=True)
@click.option("--outfileprefix", "-p",
              help="prefix to the output files",
              type=click.Path(resolve_path=True),
              default=os.path.realpath("./intervals"),
              show_default=True)
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
def extract_intervals(dogenes: bool, dotrs: bool, outfileprefix: str, gtffile: str) -> None:
    """Transform a genome annotation .gtf file into a intervals .txt file required to run velocyto.

    GTF_FILE: input file

    """
    dolog = False

    # Parse arguments
    infile = gtffile
    if outfileprefix is None:
        outfileprefix = os.path.splitext(infile)[0]
    
    # Initialize the containers
    genes = {}  # type: Dict[str, vcy.Gene]
    trs = {}  # type: Dict[str, vcy.Transcript]
    trs_by_chromstrand = defaultdict(list)  # type: Dict[str, List]

    def fix_othertrset(tr: vcy.Transcript, otrset: List) -> Tuple[int, List]:
        if len(otrset) == 0:
            return 0, []
        otrnames = []
        othergene = 1  # NOTE this is mostly for transcript model analysis (should almost never happen if gene model analysis is used)
        for otr in otrset:
            if otr.get_genename() != tr.get_genename():  # NOTE IMPORTANT we should choose if it is get_genename or get_geneid
                othergene = 2
            otrnames.append(otr.get_trname())
        otrnames.sort()
        return othergene, otrnames

    headerlines: List[Any] = []

    # Precompile regex for speed
    regex_trid = re.compile('transcript_id "([^"]+)"')
    regex_trname = re.compile('transcript_name "([^"]+)"')
    regex_geneid = re.compile('gene_id "([^"]+)"')
    regex_genename = re.compile('gene_name "([^"]+)"')
    # regex_exonno = re.compile('exon_number "(([^"]+)|([^"]+))"')

    # Loop trough the .gtf file and get transcript models and assemble gene models by fusing exons
    # NOTE: the parsing could be done using pysam tabix
    logging.debug(f"Reading {infile}")
    for nth_line, line in enumerate(open(infile)):
        if line.startswith('#'):
            headerlines.append(line)
            continue
        
        fields = line.rstrip().split('\t')
        chrom, feature_class, feature_type, start, end, junk, strand, junk, tags = fields
        # Deal with possible incongruences between .gtf and bam file chromosome naming
        if "chr" in chrom[:4]:
            chrom = chrom[3:].split(".")[0]
        else:
            chrom = chrom.split(".")[0]
        if feature_type in ("exon"):
            trid = regex_trid.search(tags).group(1)
            trname = regex_trname.search(tags).group(1)
            geneid = regex_geneid.search(tags).group(1)
            genename = regex_genename.search(tags).group(1)
            # exonno = regex_exonno.search(tags).group(1)
            start = int(start)
            end = int(end)
            chromstrand = chrom + strand
            # If encountering this transcript for the first time, create the object and add it to containers
            if trname not in trs:
                tr_tmp = vcy.Transcript(trname, trid, genename, geneid, chrom, strand)  # type: vcy.Transcript
                trs[trname] = tr_tmp
                trs_by_chromstrand[chromstrand].append(tr_tmp)
            # Add the exon to the transcript
            trs[trname].add_exon(start, end)
            # If encountering the gene of this transcript for the first time, create the object and add it to containers
            if geneid not in genes:
                gene = vcy.Transcript(genename, geneid, genename, geneid, chrom, strand)
                genes[geneid] = gene
            genes[geneid].add_exon(start, end, fuse=True)  # In this case fuse the exons if overlapping

        if nth_line % 100000 == 0:
            logging.debug(f"Read {nth_line} lines of {infile}")

    logging.debug(f"Completed reading data for {len(trs)} transcripts and {len(genes)} genes")

    #  Make containers of Transcript by chromosome+strand and bins of 100KB
    binnedtrs_by_chromstrand = {}  # type: Dict[str, Dict[int, set]]
    binsize: int = vcy.BINSIZE_BP  # The bin is used to look up for other transcripts
    for chromstrand, trlist in trs_by_chromstrand.items():
        binnedtrs = defaultdict(set)  # type: Dict[int, set]
        for t in trlist:
            for trbin in range(t.get_start() // binsize, 1 + t.get_end() // binsize):
                binnedtrs[trbin].add(t)
        binnedtrs_by_chromstrand[chromstrand] = binnedtrs

    features = []  # type: List
    if dogenes:
        features.append(('gene', genes))
    if dotrs:
        features.append(('tr', trs))
    for name, modeldict in features:
        outfile = f"{outfileprefix}_{name}_ivls.txt"
        fd = open(outfile, 'w')
        fd.write(f"# Generated from {infile}\n")
        fd.write("".join(headerlines))
        fd.write("# E: interval is exon in the models, I: interval is intron in the model\n")
        fd.write("# last two fields show other transcripts that are exon and introns within the same interval\n")
        fd.write("#TrName\tTrID\tGeneName\tGeneID\tChrom\tStrand\tIntervals\n")
        models = list(modeldict.values())  # type: List[vcy.Transcript]
        models.sort()

        logging.debug("Processing %s %ss..." % (len(models), name))
        for n_tr, tr in enumerate(models):
            chromstrand = tr.get_chromstrand()
            trstart, trend = tr.get_start(), tr.get_end()
            trname, trid, trgenename, trgeneid = tr.get_trname(), tr.get_trid(), tr.get_genename(), tr.get_geneid()

            # Retrieve transcripts within the same bin
            othertrs = set()
            for trbin in range(trstart // binsize, 1 + trend // binsize):
                othertrs.update(binnedtrs_by_chromstrand[chromstrand][trbin])
            if dolog and tr.get_trname() == "Lama4-003":
                logging.debug(tr)

            # Not used right now
            # Extend the 5' or 3' end respecting the presence of other transcripts
            for now5end, extlen in [(True, vcy.EXTENSION5_LEN), (False, vcy.EXTENSION3_LEN)]:
                if extlen > 0:
                    extension_ok = True
                    extstart, extend = (trstart - extlen, trstart) if tr.is_fw() == now5end else (trend, trend + extlen)
                    for otr in othertrs:
                        if otr.get_geneid() != trgeneid and otr.get_start() < extend and otr.get_end() > extstart:
                            for oivl in otr.sorted_exons():
                                if oivl[0] < extend and oivl[1] > extstart:
                                    logging.debug(f"Extension of {trname} hindered by {otr.get_trname()}")
                                    extension_ok = False
                        if not extension_ok:
                            break
                    if extension_ok:
                        if tr.is_fw() == now5end:
                            trstart = extstart
                        else:
                            trend = extend

            # Mark up every base in transcript locus that is own or other transcript's exon/intron
            
            # Initialize a list of list (with a list for every base)
            trlen = trend - trstart
            profile = []  # type: List[List[Any]]
            # every entry of profile constains [type_own_interval, set_other_transcripts_exon_overlap, set_other_transcripts
            for i in range(trlen):
                profile.append([0, set(), set()])  # NOTE: maybe we can implement this with numpy arrays dtype=object
            
            # Loop through the intervals (intron and exons) and mark every base with its interval type
            for ivl in tr.sorted_ivls():
                ivlstart, ivlend, ivltype = ivl  # type: Tuple[int, int, int]
                if dolog and trname == "Lama4-003":  # NOTE: I might want to get rid of this dolog checks completelly
                    logging.debug("Marking", ivlstart, ivlend, ivltype)
                for i in range(ivlstart - trstart, min(trlen, ivlend - trstart)):
                    profile[i][0] = ivltype

            # check other transctipts in the window/bin
            for otr in othertrs:
                if otr == tr:
                    continue
                # if there is some possible overlap between the two transcripts
                if otr.spans_over((trstart, trend)):
                    if dolog and trname == "Lama4-003":
                        logging.debug("Overloading %s" % otr.get_trname())
                    
                    # Visit in order all the intervals and mark the overlappings
                    for oivl in otr.sorted_ivls():
                        ostart, oend, otype = oivl
                        if ostart > trend:
                            break  # this interval is downsream the tr: stop looping
                        if oend < trstart:
                            continue  # this interval is upstream the tr: not interesting
                        pidx = 1 if otype == vcy.EXON else 2
                        if dolog and trname == "Lama4-003":
                            logging.debug("              ivl:", otr.get_trname(), ostart, oend, pidx)
                        # Add the other transcript to the profile
                        for i in range(max(0, ostart - trstart), min(trlen, oend - trstart)):
                            profile[i][pidx].add(otr)
        
            # Now delineate intervals on transcript locus that are certainly own exon, certainly own intron, maybe own exon or other
            currp = [None, None, None]  # type: List[Any]
            outivls = []  # type: List[List[Any]]
            pos = 0
            currivlstart: int
            otherexons: Any
            otherintrons: Any
            ivltype_str: str
            for p in profile:
                if p != currp:
                    if currp[0] is not None:  # first time skip
                        outivls.append([currivlstart, pos + trstart, ivltype_str, otherexons, otherintrons])
                    currp = p
                    ivltype_str = 'E' if currp[0] == vcy.EXON else 'I'
                    otherexonstatus, otherexons = fix_othertrset(tr, p[1])
                    ivltype_str += '-eE'[otherexonstatus]
                    otherintronstatus, otherintrons = fix_othertrset(tr, p[2])
                    ivltype_str += '-iI'[otherintronstatus]
                    currivlstart = pos + trstart
                pos += 1
            if currp[0] is not None:  # do it a last time for the last iteration
                outivls.append([currivlstart, pos + trstart, ivltype_str, otherexons, otherintrons])
        
            # Write output file
            fd.write("%s\t%s\t%s\t%s\t%s\t%s" % (trname, trid, trgenename, trgeneid, tr.get_chrom(), tr.get_strand()))
            for ivl in outivls:
                fd.write("\t%s-%s:%s:%s:%s" % (ivl[0], ivl[1], ivl[2], ",".join(ivl[3]), ",".join(ivl[4])))
            fd.write("\n")

            if (n_tr + 1) % 2500 == 0:
                logging.debug(f"{n_tr+1} {name}s already processed")

        fd.close()

    outputs_names = ""
    for name, modeldict in features:
        outputs_names += f"{outfileprefix}_{name}_ivls.txt"
    logging.debug(f"Done without errors. Find the outputs at {outputs_names}")
