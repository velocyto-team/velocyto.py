.. _changelog:

=========
Changelog
=========

* :release:`0.17.8 <2018-05-30>`
* :bug:`-` Avoids an error with an older versions of cellranger not outputting the _log file
* :release:`0.17.6 <2018-05-03>`
* :feature:`-` Be more forgiving on the requirements of the gtf file: `gene_name`, `transcipt_name` and `exon_number` are not required fields anymore
* :support:`-` Add section of the docs, explaining better what the bam and gtf files are expected to contain
* :release:`0.17.5 <2018-04-19>`
* :feature:`-` Add `filter_genes_by_phase_portrait` as the main api to filter genes supporting versatile set of conditions
* :release:`0.17.2 <2018-04-08>`
* :bug:`-` Fix an error when aligners report consecutive insertion/deletion instead of mismatch 
* :release:`0.17.0 <2018-04-08>`
* :feature:`-` Add subcommand ``run_dropest`` as a shortcut to run dropEst preprocessed data (including Indrops and DropSeq) 
* :support:`-` Deprecation warning on `default` functions: they were being misused by the users.
* :feature:`-` Add the set of subcommands ``veloctyo tools`` to bridge velocyto with other software (for now DropEst)
* :feature:`-` Add cosine projection penalty
* :feature:`-` Change behavior in no-barcode-list mode: use a very permissive heuristic of < 80 molecules per cell as the threshold to thrash empty droplets / no-cell events
* :bug:`-` Fix a skip repeat error with SmartSeq2 pipeline
* :feature:`-` run automatically randomized (negative) control for the velocity. Added plotting options for the randomized control visualization
* :bug:`-` Fix colormap bug with matplotlib 2.2.0
* :feature:`-` Add optional feature selection for the unspliced
* :feature:`-` Add size factor normalization option
* :feature:`-` Add possibility to constraint knn averaging: when turned on avoids edges between cells of specified groups
* :bug:`-` Fix an error in filter_cells: colors array is now filtered as well
* :feature:`-` Improve the debug molecular report option to support hdf5
* :feature:`-` It supports SmartSeq2 and has a new run_smartseq2 command
* :feature:`-` It supports multiple bam files input, it can interpret the file(s) as either as one-file-one-cell or just as batches to be analyzed together. IMPORTANT: a cell cannot be distributed over different bamfiles!
* :feature:`-` Generalized logic to include more layers than just Spliced, Unspliced, Ambiguous
* :feature:`-` ``--without-umi`` option allows analyzing UMI-less data such as SmartSeq2
* :feature:`-` support different verbosity levels with the ``-v`` flag
* :support:`-` ``--multimap`` option was removed because it could have yield incorrect results depending on the output format chosen for the aligner
* :feature:`-` Support barcodes of different length
* :feature:`-` Improved compatibility with InDrops and STRT through the ``--umi-extentions``. It allows the same pipeline to be applied to methods with short molecular barcode that cannot be used call a unique molecule without the gene mapping information.
* :release:`0.13.5 <2018-02-07>`
* :bug:`-` Fix a bug that caused extremely slow runtimes when the input bam was not position sorted. Now `velocyto` will raise an error and ask the user to sort the file using samtools.
* :support:`-` Improve the changelog structure
* :release:`0.13.4 <2018-01-25>`
* :bug:`-` A change in slicing related to an API change of `__getattr__` in loompy2 
* :release:`0.13.3 <2018-01-25>`
* :bug:`-` Catch another error due to the API change of `.create` in loompy2 
* :release:`0.13.2 <2018-01-25>`
* :bug:`-` Catch error due to the API change of `.create` in loompy2 
* :bug:`-` Fix an incompatibility with loompy2 related to column and row attributes changing from dict to an object
* :release:`0.13.1 <2018-01-22>`
* :feature:`-` Sample metadata file can be specified with different csv formats (the format will be determined automatically)
* :release:`0.13.0 <2018-01-19>`
* :bug:`-` Sometimes velocyto missed to detect and warn the user that the `.gtf` genome annotation file was not sorted, this could have caused undetected errors in the analysis. If you run velocyto without sorting the .gtf, we suggest rerunning.
* :feature:`-` CLI does not require presorting the gtf files. To reduce possibility of incorrect usage, now .gtf file sorting sorting is performed in memory (and not saved).
* :feature:`-` Improve documentation: remove information about sorting .gtf files. This procedure is not needed anymore.
* :release:`0.12.4 <2018-01-18>`
* :bug:`40` Error in hdf5 serialization when using cluster label as object array is now fixed
* :release:`0.12.3 <2018-01-17>`
* :bug:`38` Error in running run10x is now fixed
* :release:`0.12.2 <2018-01-12>`
* :bug:`37` Initial cell size array gets updated properly when filtering cells
* :release:`0.12.1 <2018-01-04>`
* :bug:`35` Attempted to fix a reported bug running velocyto CLI
* :release:`0.12.0 <2017-12-17>`
* :feature:`-` Add possibility to export pickle containing information of every molecule
* :feature:`-` Remove the subcommand ``multi10x``
* :bug:`- major` Incorrect 0-based indexing for splicing junction corrected (was not causing problems because buffered by MIN_FLANK) 
* :bug:`- major` Many small bug fixes
* :bug:`31 major` Memory usage bug should be solved.
* :feature:`-` Large parts of the documentation rewritten to match the changes in API
* :feature:`-` Make the CLI simpler removing the extract interval step. 
  Now the source .gtf files can be provided directly, they should be provided sorted using :code:`sort -k1,1 -k7,7 -k4,4n -o [OUTFILE] [INFILE]`
* :feature:`-` Changelog added to the doc
* :support:`-` Update the documentation for the new  :ref:`CLI <cli>`
* :feature:`-` Support different Logic levels
* :feature:`-` Pipeline now consider all the possible transcript models that could be supported by a set of reads individually and then decides on the spliced/unspliced/ambiguous count.
* :release:`0.11.0 <2017-12-01>`
* :bug:`- major` fix a bug with ambiguous molecules counting and version bump
* :release:`0.10.3 <2017-11-23>`
* :bug:`- major` The debug and sampleid option had the same short flag `-d`
* :release:`0.10.2 <2017-11-18>`
* :release:`0.10.1 <2017-11-18>`
* :feature:`-` further ~5x speedup of cython functions making them 100% C and using malloc instead of memory views
* :release:`0.10.0 <2017-11-18>`
* :feature:`-` Add support for DropSeq pipelines where the barcode flags in the bam file are `XC` and `XM` instead of `CB` and `UB`
* :bug:`- major` Using sphinx 1.7 sorts the autodoc API correctly
* :release:`0.9.13 <2017-11-04>`
* :release:`0.9.12 <2017-11-04>`
* :release:`0.9.11 <2017-11-03>`
* :feature:`-` Improve the docs
* :release:`0.9.10 <2017-11-02>`
* :release:`0.9.9 <2017-10-31>`
* :release:`0.9.8 <2017-10-26>`
* :release:`0.9.7 <2017-10-25>`
* :release:`0.9.6 <2017-10-25>`