.. _cli:

CLI Usage Guide
===============

Introduction
------------

After you have velocyto correctly installed on your machine (see :ref:`installation tutorial <install>`) the ``velocyto`` command will become available in the terminal.
``velocyto`` is a command line tool with subcomands. You can get info on all the available commands typing ``velocyto --help``. You will get the following output:

.. include:: ../substitutions/velocyto.txt


You can further query for information on each subcommand by typing ``velocyto COMMANDNAME --help``.

Alternatively you can visit the online :ref:`api description page <cliapi>` that includes usage information for all the subcommands.

Preparation
-----------

Download genome annotation file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download a genome annotation (.gtf file) for example from `GENCODE <http://www.gencodegenes.org/>`_ or `Ensembl <http://www.ensembl.org/info/data/ftp/index.html>`_. If you use the  ``cellranger`` pipeline, you should download the gtf that comes prepackaged with it `here <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references>`_.

Download expressed repeats annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. note::
   This step is optional.

You might want to mask expressed repetitive elements, since those count could constitute a confounding factor in the downstream analysis.
To do so you would need to download an appropriate expressed repeat annotation (for example `from UCSC genome browser <https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf>`_ and **make sure to select GTF**  as output format).


Running ``velocyto``
--------------------

The general purpose command to start the pipeline for read counting is ``velocyto run``.
The ``run`` defaults are appropriate for the analysis of both 10X Genomics v1/v2 and InDrops 3' chemistry.

A typical use of ``run`` is:

::

    velocyto run -b filtered_barcodes.tsv -o output_path -m repeat_msk_srt.gtf possorted_genome_bam.bam mm10_annotation.gtf


The general signature for the ``run`` subcommand is:

.. include:: ../substitutions/run.txt

.. note::
    Execution time is ~3h for a typical sample but might vary significantly by sequencing depth and cpu power. 
    
.. warning::
   Running velocyto without specifying a filtered barcode set (:code:`-b/--bcfile` option) is not recommended, do it at your own risk. In thsi way, the counter will use all the cell barcodes it encounters. It might result in long runtimes, large memory alloactions and big output matrix.


Notes on ``velocyto run``
~~~~~~~~~~~~~~~~~~~~~~~~~

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

- Permissive10X
- ValidatedIntrons10X (*Default)
- Stricter10X
- ObservedSpanning10X

Despite the name (that designates their original design for the 10X platform) the logics generalize well to similar chemistries (e.g. Drop-seq).

.. hint::
    Custom logics supporting peculiarities of other chemistries can be implemented simply by creating a class that inherits from :ref:`Logic <logicapi>`.

Run on a single or multiple 10X Chromium samples
------------------------------------------------

``velocyto`` supports a shortcut to run directly on one or more `cellranger` output folders (e.g. this is the folder containing the subfolder: ``outs``, ``outs/analys`` and ``outs/filtered_gene_bc_matrices``).

For example if we want to run the pipeline on the folder ``mypath/sample01``. We would do:

::

    velocyto run10x -m repeat_msk_srt.gtf mypath/sample01 mm10_annotation.gtf

The full signature of the command is:

.. include:: ../substitutions/run10x.txt


About the output .loom file
---------------------------

The main result file is a 4-layered `loom file <http://loompy.org/loompy-docs/format/index.html>`_ : `sample_id.loom`. 

A valid .loom file is simply an HDF5 file that contains specific groups representing the main matrix as well as row and column attributes.
Because of this, .loom files can be created and read by any language that supports HDF5. 

.loom files can be easily handled using the `loompy package <http://loompy.org>`_.

Get started with the analysis
-----------------------------

At this point you are ready to start analyzing your ``.loom`` file. To get started read our :ref:`analysis tutorial <analysis>` and have a look at the :ref:`notebooks examples <notebooks>`.
