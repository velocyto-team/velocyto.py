
from velocyto.commands._run import _run
# from velocyto.commands.run10x import run10x
from velocyto.commands.common import logicType, init_logger
from pathlib import Path
from loguru import logger
import numpy as np
import sys
sys.path[0] = str(Path(sys.path[0]).parent)
# from memory_profiler import profile, LogFile


# run10x(
#     samplefolder=Path("/mnt/vault/PAM/PAM1"),
#     gtffile=Path("/mnt/vault/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz"),
#     mask=Path("/mnt/group/references/genomic/homo_sapiens/sequences/grch3810_repeat_mask.gtf")
# )

def main():
    logger.remove()
    init_logger(3, msg_format="<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>")
    # sys.stdout = LogFile('memory_profile_log', reportIncrementFlag=False)

    samplefolder = Path('/mnt/vault/workspace/analysis/pam')

    umap_file = samplefolder.joinpath("outs/per_sample_outs/analysis/umap/gene_expression_2_components/projection.csv")
    umap = np.loadtxt(umap_file, usecols=(1, 2), delimiter=",", skiprows=1)

    clusters_file = samplefolder.joinpath("outs/per_sample_outs/analysis/clustering/gene_expression_graphclust/clusters.csv")
    labels = np.loadtxt(clusters_file, usecols=(1,), delimiter=",", skiprows=1)


    additional_ca = {
        "_X": umap[:, 0].astype("float32"),
        "_Y": umap[:, 1].astype("float32"),
        "Clusters": labels.astype("int") - 1
    }


    _run(
            bamfile=samplefolder.joinpath('subsample_alignments.bam'),
            gtffile=samplefolder.joinpath('genes.gtf'),
            bcfile=Path("/mnt/group/software/cellranger-7.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt"),
            outputfolder=samplefolder.joinpath('velocyto'),
            sampleid="PAM1",
            metadatatable=None,
            repmask=samplefolder.joinpath('grch3810_repeat_mask.gtf'),
            onefilepercell=False,
            logic=logicType.Permissive10X,
            without_umi=False,
            umi_extension="no",
            multimap=False,
            test=False,
            samtools_threads=4,
            samtools_memory="4G",
            dump="0",
            loom_numeric_dtype='uint32',
            verbose=True,
            bughunting=True,
            additional_ca=additional_ca,
            is_10X=True,
            samplefolder=samplefolder,
        )

main()
