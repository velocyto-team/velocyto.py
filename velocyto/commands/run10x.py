import sys
import os
import glob
import re
import gzip
import click
import array
import loompy
import numpy as np
import csv
from collections import defaultdict
import logging
from typing import *
import velocyto as vcy
from ._run import _run


logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


@click.command(short_help="Runs the velocity analysis for a Chromium Sample")
@click.argument("samplefolder",
                type=click.Path(exists=True,
                                file_okay=False,
                                dir_okay=True,
                                readable=True,
                                writable=True,
                                resolve_path=True))
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.option("--metadatatable", "-s",
              help="Table containing metadata of the various samples (csv fortmated rows are samples and cols are entries)",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--repmask", "-m",
              help=".gtf file containing intervals to mask",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--logic", "-l",
              help="The logic to use for the filtering (default: Default)",
              default="Default")
@click.option("--multimap", "-M",
              help="Use reads that did not map uniquely (default: False)",
              default=False,
              is_flag=True)
@click.option("--samtools-threads", "-@",
              help="The number of threads to use to sort the bam by cellID file using samtools",
              default=16)
@click.option("--samtools-memory",
              help="The number of MB used for every thread by samtools to sort the bam file",
              default=2048)
def run10x(samplefolder: str, gtffile: str,
           metadatatable: str, repmask: str, logic: str, samtools_threads: int, samtools_memory: int,
           multimap: bool) -> None:
    """Runs the velocity analysis for a Chromium 10X Sample

    10XSAMPLEFOLDER specifies the cellranger sample folder

    GTFFILE genome annotation file
    """

    # Check that the 10X analysis was run successfully
    if "Pipestance completed successfully!" not in open(os.path.join(samplefolder, "_log")).read():
        raise IOError("The outputs are not ready")
    bamfile = os.path.join(samplefolder, "outs", "possorted_genome_bam.bam")
    bcfile = glob.glob(os.path.join(samplefolder,
                       os.path.normcase("outs/filtered_gene_bc_matrices/*/barcodes.tsv")))[0]
    outputfolder = os.path.join(samplefolder, "velocyto")
    sampleid = os.path.basename(samplefolder.rstrip("/").rstrip("\\"))
    assert not os.path.exists(os.path.join(outputfolder, f"{sampleid}.loom")), "The output already exist. Aborted!"
    additional_ca = {}
    try:
        tsne_file = os.path.join(samplefolder, "outs", "analysis", "tsne", "2_components", "projection.csv")
        if os.path.exists(tsne_file):
            tsne = np.loadtxt(tsne_file, usecols=(1, 2), delimiter=',', skiprows=1)
            additional_ca["_X"] = tsne[:, 0].astype('float32')
            additional_ca["_Y"] = tsne[:, 1].astype('float32')

        clusters_file = os.path.join(samplefolder, "outs", "analysis", "clustering", "graphclust", "clusters.csv")
        if os.path.exists(clusters_file):
            labels = np.loadtxt(clusters_file, usecols=(1, ), delimiter=',', skiprows=1)
            additional_ca["Clusters"] = labels.astype('int') - 1

    except Exception:
        logging.error("Some IO problem in loading cellranger tsne/pca/kmeans files occurred!")

    return _run(bamfile=bamfile, gtffile=gtffile, bcfile=bcfile, outputfolder=outputfolder,
                sampleid=sampleid, metadatatable=metadatatable, repmask=repmask,
                logic=logic, molrep=False, multimap=multimap, test=False, samtools_threads=samtools_threads,
                samtools_memory=samtools_memory, additional_ca=additional_ca)
