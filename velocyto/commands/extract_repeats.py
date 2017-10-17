import logging
import copy
import os
import re
import sys
import click
import gzip
import bisect
import subprocess
from collections import defaultdict
from typing import *
import velocyto as vcy

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


@click.command(short_help="Transform a repeats .gtf file into a intervals .gtf file.")
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.option("--tolerance", "-t",
              default=5,
              show_default=True)
def extract_repeats(gtffile: str, tolerance: int) -> None:
    """Transform a repeats .gtf file into a intervals .gtf file.

    GTFFILE  .gtf file with annotated repeats to mask
    """
    
    sorted_filename = gtffile.split(".")[0] + "_sorted.gtf"
    logging.debug(f"Sorting by `sort -k1,1 -k7,7 -k4,4n {gtffile} > {sorted_filename}`")
    with open(sorted_filename, "w") as f:
        p1 = subprocess.run(["sort", "-k1,1", "-k7,7", "-k4,4n", gtffile],
                            stdout=f)
        
    # Parse arguments
    logging.debug(f'Opening {sorted_filename}, output will be at {gtffile.split(".")[0] + "_joined.gtf"}')
    fin = open(sorted_filename)
    fout = open(gtffile.split(".")[0] + "_joined.gtf", "w")

    curr_chromstrand = ""
    curr_end = -1000
    TOLERANCE = tolerance
    headerlines = []
    for line in fin:
        if line.startswith('#'):
            headerlines.append(line)
            continue
        
        fields = line.rstrip().split('\t')
        chrom, feature_class, feature_type, start, end, junk, strand, junk, tags = fields
        start = int(start)
        end = int(end)
        chromstrand = chrom + strand
        
        curr_chrom: str
        curr_feature_class: str
        curr_feature_type: str
        curr_n: int
        curr_start: int
        curr_tags: str
        # If there is no overlap with the next
        if chromstrand != curr_chromstrand or start > curr_end + TOLERANCE:
            if curr_chromstrand != "":
                new_entry = (curr_chrom, curr_feature_class, curr_feature_type,
                             str(curr_start), str(curr_end), ".", strand, ".",
                             f"number_elements {curr_n}; {curr_tags}")
                fout.write("\t".join(new_entry).strip() + "\n")
            
            curr_chrom = chrom
            curr_feature_class = feature_class
            curr_feature_type = feature_type
            curr_start = start
            curr_end = end
            curr_n = 1
            curr_strand: str = strand
            curr_tags = tags
            curr_chromstrand = chromstrand
        else:
            gap = start - curr_end
            curr_end = end
            curr_n += 1
            curr_tags = f"{curr_tags} gap {gap}; {tags}" if gap > 0 else curr_tags + tags
    fout.close()
    fin.close()

    logging.debug(f"Deleting the sorted file {sorted_filename}")
    os.remove(sorted_filename)
