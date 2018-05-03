.. _cli:

CLI Usage Guide
===============

Introduction
------------

After you have velocyto correctly installed on your machine (see :ref:`installation tutorial <install>`) the ``velocyto`` command will become available in the terminal.
``velocyto`` is a command line tool with subcomands. You can get quick info on all the available commands typing ``velocyto --help``. You will get the following output:

.. include:: ../substitutions/velocyto.txt


More detailed information can be obtained using ``--help`` for each subcommand (by typing ``velocyto COMMANDNAME --help``).

Alternatively you can visit the online :ref:`api description page <cliapi>` that includes usage information for all the subcommands.

Preparation
-----------

Download genome annotation file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download a genome annotation (.gtf file) for example from `GENCODE <http://www.gencodegenes.org/>`_ or `Ensembl <http://www.ensembl.org/info/data/ftp/index.html>`_. If you use the  ``cellranger`` pipeline, you should download the gtf that comes prepackaged with it `here <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references>`_.

Download expressed repeats annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note::
   This step is optional.

You might want to mask expressed repetitive elements, since those count could constitute a confounding factor in the downstream analysis.
To do so you would need to download an appropriate expressed repeat annotation (for example `from UCSC genome browser <https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf>`_ and **make sure to select GTF**  as output format).


Running ``velocyto``
--------------------

The general purpose command to run the read counting pipeline is ``velocyto run``.
However, for some of the most commonly used scRNA-seq chemistries, we provide a set of ready-to-use subcommands.
The currently available are: ``run10x``, ``run_smartseq2``, ``run_dropest``
These subcommands are just wrappers of the main command ``velocyto run``. They take care of passing the appropriate options for each technology; furthermore they performa a minimal check the inputs provided make sense and some of them infer the path of some of the input files.

Please regard these subcommands as the recommended and easiest way to run `velocyto`, especially if you are unsure of all the options to pass to ``velocyto run``.
For more flexibility and advanced usage, ``velocyto run`` should be used directly.
Furthermore, to adapt velocyto to custom/new techniques the user may want to consider modification of the counting pipeline, this does not require deep rewrite of the internals but just the creation of a new `logic`, for more information consult the section about the :ref:`Logic interface API <logicapi>`.

We will now describe the use of the technique specific subcommands.
If you are interested in running velocyto with only one technique you can directly jump to that section without loss of information (this also means that some of the information will be repeated).

``run10x`` - Run on 10X Chromium samples
----------------------------------------

``velocyto`` includes a shortcut to run the counting directly on one or more `cellranger` output folders (e.g. this is the folder containing the subfolder: ``outs``, ``outs/analys`` and ``outs/filtered_gene_bc_matrices``).

The full signature of the command is:

.. include:: ../substitutions/run10x.txt

For example if we want to run the pipeline on the `cellranger` output folder ``mypath/sample01``. We would do:

::

    velocyto run10x -m repeat_msk.gtf mypath/sample01 somepath/refdata-cellranger-mm10-1.2.0/genes/genes.gtf

Where ``genes.gtf`` is the genome annotation file provided with the cellranger pipeline.
``repeat_msk.gtf`` is the repeat masker file described in the `Preparation` section above.

.. note::
    Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power. 


``run_smartseq2`` - Run on SmartSeq2 samples
--------------------------------------------

``velocyto`` includes a shortcut to perform the read counting for UMI-less, not stranded, full-length techniques such as SmartSeq2.

The full signature of the command is:

.. include:: ../substitutions/run_smartseq2.txt

Typically SmartSeq2 bam files are generated and organized by well/cell in a folder structure similar to the following

::

    plateX/A01/A01.bam
    plateX/A02/A02.bam
    plateX/A03/A03.bam
    ...

For this reason, ``run_smartseq2`` command accepts multiple inputs so that all the files can be analyzed in one run and all output stored in a single `.loom` file. 
For example if we want to run the pipeline on a SmartSeq2 plate whose folder structure is organized as above we would just use wild-card (glob) expansion as follows.

::

    velocyto run_smartseq2 -o OUTPUT -m repeat_msk.gtf -e MyTissue plateX/*/*.bam mm10_annotation.gtf

Where ``mm10_annotation.gtf`` and ``repeat_msk.gtf`` are the genome annotation and repeat masker files described in the `Preparation` section above.
Finally note that the output `.loom` file in that case will have an extra layer "spanning".

.. note::
    The input bam files need to be sorted by position, if this is not already the case, this is simply achieved running `samtools sort A01.bam -o sorted_A01.bam`. 

.. note::
    Execution time might vary significantly by sequencing depth and cpu power but usually does not exceed the 6h for a typical sample 


``run_dropest`` - Run on DropSeq, InDrops and other techniques
--------------------------------------------------------------

``velocyto`` includes a shortcut to perform molecule counting on all the techniques supported by the DropEst pipeline, this includes different versions of DropSeq and InDrops.
This is particularly convenient since the output from the pipeline is similar for different techniques allowing the use of a single command.

.. note::
    If you prefer using another pipeline or technique not supported by DropEst, note that you can still use the core command ``velocyto run`` but no shortcut is provided yet.
    We are eager to work with implement more shortcut for other pipelines and techniques (please contact us if you are the developer or can help us integrate velocyto seamlessly in other pipelines).

The full signature of the command is:

.. include:: ../substitutions/run_dropest.txt

Let's start with a full example including the steps to run DropEst correctly.
Explaining how to run DropEst is outside the scope of this documentation. For more information on installation and usage of DropEst refer to its documentation.
I will assume I am analyzing an InDrops sample downloaded from SRA.

First of all set 10 as minimum barcode quality. Then I start by calling the `droptag` command as follows:

::

    ./droptag -c ./configs/indrop_v1_2.xml ~/mydata/SRR5945694_2.fastq.gz ~/mydata/SRR5945694_1.fastq.gz

The 

::

    STAR --genomeDir ~/cellranger/refdata-cellranger-mm10-1.2.0/star/ \
           --readFilesIn SRR5945695_1.fastq.gz.tagged.1.fastq.gz,SRR5945695_1.fastq.gz.tagged.2.fastq.gz,SRR5945695_1.fastq.gz.tagged.3.fastq.gz,SRR5945695_1.fastq.gz.tagged.4.fastq.gz,SRR5945695_1.fastq.gz.tagged.5.fastq.gz,SRR5945695_1.fastq.gz.tagged.6.fastq.gz,SRR5945695_1.fastq.gz.tagged.7.fastq.gz,SRR5945695_1.fastq.gz.tagged.8.fastq.gz \
           --outSAMmultNmax 1 \
           --runThreadN 6 \
           --readNameSeparator space \
           --outSAMunmapped Within \
           --outSAMtype BAM SortedByCoordinate \
           --outFileNamePrefix SRR5945695_1 \
           --readFilesCommand gunzip -c

Note the -m option will read the config_desc.xml file that should have the barcodes_file option correctly selected.
For InDrops for example this would be as following:

.. code-block:: xml
    
    <barcodes_file>[YOURPATH]/dropEst/data/barcodes/indrop_v1_2</barcodes_file>

And then I run DropEst core command:

::

    dropest -m -V -b \
              -o ~/mydata/SRR5945695/SRR5945695_dropEst \
              -g ~/cellranger/refdata-cellranger-mm10-1.2.0/genes/genes.gtf \
              -L eiEIBA \
              -c ~/mysource/dropEst/configs/config_desc.xml \
              ~/mydata/SRR5945695_1/SRR5945695_1Aligned.sortedByCoord.out.bam

Then I run ``velocyto tools dropest_bc_correct`` with the bamfile as the only argument. This will:
(1) Make a new bam with the barcodes substituted with the corrected ones, taking this info from the dropEst R dump
(2) Generate the required file containing the allowed barcodes

The bam file outputed by dropEst does not contain error-corrected but raw cell barcodes so we will have to make a new corrected bam file using the infromation otuputed.
(Future version of DropEst will output the error-corrected barcodes).

To do that run:

::

    velocyto tools dropest_bc_correct ~/mydata/SRR5945695_1/SRR5945695_1Aligned.sortedByCoord.out.tagged.bam

Finally I run velocyto:

::

    velocyto run_dropest -o ~/mydata/SRR5945695_results -m rep_mask.gtf ~/mydata/SRR5945695_1/correct_SRR5945695_1Aligned.sortedByCoord.out.tagged.bam mm10_annotation.gtf


``run`` - Run on any technique (Advanced use)
---------------------------------------------

The general signature for the ``run`` subcommand is:

.. include:: ../substitutions/run.txt


A typical use of ``run`` is:

::

    velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf possorted_genome_bam.bam mm10_annotation.gtf



.. note::
    The input bam file needs to be sorted by position, this can be achieved running `samtools sort mybam.bam -o sorted_bam.bam`. In cellranger generated bamfiles are already sorted this way. 

.. note::
    Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power. 
    
.. warning::
   Running velocyto without specifying a filtered barcode set (:code:`-b/--bcfile` option) is not recommended, do it at your own risk. In this way, the counter will use all the cell barcodes it encounters. It might result in long runtimes, large memory allocations and big output matrix.


Notes on first runtime and parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As one of its first steps ``velocyto run`` will try to create a copy of the input .bam files sorted by cell-barcode. The sorted .bam file will be placed in the same directory as the original file and it will be named ``cellsorted_[ORIGINALBAMNAME]``.
The sorting procedure uses ``samtools sort`` and it is expected to be time consumning, because of this, the procedurre is perfomed in parellel by default. It is possible to control this parallelization using the parameters ``--samtools-threads`` and  ``--samtools-memory``.

.. note::
    If the file ``cellsorted_[ORIGINALBAMNAME]`` exists, the sorting procedure will be skipped and the file present will be used.

.. warning::
    Most of the ``velocyto`` pipeline is single threaded and several instances can be run on the same multicore machine to process your samples in a time effective way. 
    However, because of the above mentioned multithreaded call to ``samtools sort``, running several instances of ``veloctyo run`` might end up using the memory and cpu of your system and possbily result in runtime errors.
    Therefore for batch jobs we suggest to first call :code:`samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam` sequentially and only then running ``velocyto`` 

Run with different logics
~~~~~~~~~~~~~~~~~~~~~~~~~

The rules used to call spliced, unspliced and ambiguous molecules from the reads mappings can be set using the :code:`--logic` parameter.
The behavior of the counter can be modified using one of the different logics supported.
Every logic has a different sensitivity. The currently available are:

- *Default*
- SmartSeq2
- ValidatedIntrons10X 
- Stricter10X
- ObservedSpanning10X

Despite some of the names (that is kept for backwards compatibility and designates their original design for the 10X platform) the logics are tested and work for similar chemistries (e.g. InDrops).

.. hint::
    Custom logics supporting peculiarities of other chemistries can be implemented simply by creating a class that inherits from :ref:`Logic <logicapi>`.


Requirements on the input files
-------------------------------

velocyto assumes that the ``bam`` file that is passed to the CLI contains a set of information and that some upstream analysis was performed on them already.
In particular the ``bam`` file will have to:

1. Be sorted by mapping position.
2. Represents either a single sample (multiple cells prepared using a certain barcode set in a single experiment) or single cell.
3. Contain an `error corrected` cell barcodes as a TAG named ``CB`` or ``XC``.
4. Contain an `error corrected` molecular barcodes as a TAG named ``UB`` or ``XM``.

.. note::
    For SmartSeq2 bam files (3) and (4) are not required because it consists of one bam file per cell and no umi are present.


velocyto assumes that the ``gtf`` file follows the `GENCODE gtf format description <https://www.gencodegenes.org/gencodeformat.html>`_.
Hoever some mandatory field are relaxed to extend compatibility to a wider set of gtf files.
In particular the ``gtf`` file will have to:

1. Contain the 3rd column entry ``feature-type``. Note that only the `exon` entry of the gtf file marked as `exon`in this column will be considered and therefore the requirements below only apply to the ``exon`` labeled lines.
2. Contain, in the 9th column, the key-value pair ``transcript_id``, containing an unique identified for the transcript model. 
3. Contain, in the 9th column, the key-value pair ``transcript_name`` (Optional, if not present it will be set to the value of `transcript_id`)
4. Contain, in the 9th column, the key-value pair ``gene_id``, containing an unique identified for the gene. 
5. Contain, in the 9th column, the key-value pair ``gene_name`` (Optional, if not present it will be set to the value of `gene_id`)
6. Contain, in the 9th column, the key-value pair ``exon_number`` (Reccomended but optional, if not provided velocyto will sort exons in memory and number them)

.. warning::
    Annotation of artificial chromosomes such as the ones generated to count ERCC spikes or transgenes (GFP, Tomato, etc.) need also to contain the information above.


About the output .loom file
---------------------------

The main result file is a 4-layered `loom file <http://loompy.org/loompy-docs/format/index.html>`_ : `sample_id.loom`. 

A valid .loom file is simply an HDF5 file that contains specific groups representing the main matrix as well as row and column attributes.
Because of this, .loom files can be created and read by any language that supports HDF5. 

.loom files can be easily handled using the `loompy package <http://loompy.org>`_.

Merging multiple samples/lanes in a single file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The merging of different samples/lanes in the same loom file can be performed simply using the ``loompy`` library.
This is usually just a single line:


.. code-block:: python

    loompy.combine(files, output_filename, key="Accession")

or if you want more control on the exact memory that is allocated or subset your data before merging you can do something like:

.. code-block:: python

    
    files = ["file1.loom","file2.loom","file3.loom","file4.loom"]
    # on the command line do: cp file1.loom merged.loom 
    ds = loompy.connect("merged.loom")
    for fn in files[1:]:
        ds.add_loom(fn, batch_size=1000)

Get started with the analysis
-----------------------------

At this point you are ready to start analyzing your ``.loom`` file. To get started read our :ref:`analysis tutorial <analysis>` and have a look at the :ref:`notebooks examples <notebooks>`.
