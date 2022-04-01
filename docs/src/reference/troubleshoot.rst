Troubleshoot common issues
==========================

If you have a question that is not addressed here, please don't hestitate to `ask for help <https://discussion.nextstrain.org/>`__

My country / division does not show up on the map
-------------------------------------------------

This is most often a result of the country / division not being present in `the file defining the latitude & longitude of each deme <https://github.com/nextstrain/ncov/blob/master/defaults/lat_longs.tsv>`__. Adding it to that file (and rerunning the Snakemake rules downstream of this) should fix this.

My trait (e.g. division) is grey instead of colored
---------------------------------------------------

We generate the colors from the ``colors`` rule in the Snakefile, which uses the `ordering TSV <https://github.com/nextstrain/ncov/blob/master/defaults/color_ordering.tsv>`__ to generate these. See :doc:`../guides/workflow-config-file` for more info.

*A note about locations and colors:* Unless you want to specifically override the colors generated, it's usually easier to *add* information to the default ``ncov`` files, so that you can benefit from all the information already in those files.

My genomes aren't included in the analysis
------------------------------------------

There are a few steps where sequences can be removed:

-  During the ``filter`` step:

   -  Samples that are included in `the exclude file <https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt>`__ are removed
   -  Samples that fail the current filtering criteria, as defined in the ``parameters.yaml`` file, are removed. You can modify the snakefile as desired, but currently these are:

      -  Minimum sequence length of 25kb
      -  No ambiguity in (sample collection) date

   -  Samples may be randomly removed during subsampling; see :doc:`../guides/workflow-config-file` for more info.
   -  During the ``refine`` step, where samples that deviate more than 4 interquartile ranges from the root-to-tip vs time are removed

Sequencing and alignment errors
-------------------------------

Genome sequencing, bioinformatic processing of the raw data, and alignment of the sequences are all steps were errors can slip in. Such errors can distort the phylogenetic analysis. To avoid sequences with known problems to mess up the analysis, we keep a list of problematic sequences in ``config/exclude.txt`` and filter them out. To facilitate spotting such problematic sequences, we added an additional quality control step that produces the files

-  ``results/sequence-diagnostics.tsv``
-  ``results/flagged-sequences.tsv``
-  ``results/to-exclude.txt``

These files are the output of ``scripts/diagnostics.py`` and are produced by rule ``diagnostic``. The first file contains statistics for every sequence in the alignment, sorted by divergence worst highest to lowest. The second file contains only those sequences with diagnostics exceeding thresholds each with their specific reason for flagging – these are sorted by submission date (newest to oldest). The third file contains only the names of the flagged sequences and mirrors the format of ``config/exclude.txt``. These names could be added to ``config/exclude.txt`` for permanent exclusion. Note, however, that some sequences might look problematic due to alignment issues rather than intrinsic problems with the sequence. The flagged sequences will be excluded from the current run.

To only run the sequence diagnostic, you can specify any of the three above files as target, or use the ``diagnostic`` target:

.. code:: bash

   nextstrain build ... diagnostic

In addition, we provide rules to re-examine the sequences in ``config/exclude.txt``. By running

.. code:: bash

   nextstrain build ... diagnose_excluded

the workflow will produce

-  ``results/excluded-sequence-diagnostics.tsv``
-  ``results/excluded-flagged-sequences.tsv``
-  ``results/check-exclusion.txt``

These files are meant to facilitate checking whether sequences in ``config/exclude.txt`` are excluded for valid reasons.
