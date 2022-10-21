Remote inputs
=============

This page provides an overview of intermediate files which Nextstrain produces via daily workflow runs. Where appropriate, these files can be starting points for the `ncov workflow <https://github.com/nextstrain/ncov/>`__ (discussed below).

We have two GitHub repositories which routinely upload files to `S3 buckets <https://aws.amazon.com/s3/>`__: `ncov-ingest <https://github.com/nextstrain/ncov-ingest/>`__ and `ncov <https://github.com/nextstrain/ncov/>`__. Each of those runs separate pipelines for GISAID and GenBank (aka "open") data sources; these pipelines start with data curation and QC steps and end with the phylogenetic analyses you can see on `nextstrain.org <https://nextstrain.org/sars-cov-2/>`__

The GISAID data is stored at ``s3://nextstrain-ncov-private`` and is not publicly available, in line with the GISAID Terms of Use (this is used internally by Nextstrain).

The open (GenBank) data is publicly available at three endpoints:

-  ``https://data.nextstrain.org/files/ncov/open/``
-  ``s3://nextstrain-data/files/ncov/open/``
-  ``gs://nextstrain-data/files/ncov/open/`` (mirrored daily from S3 by the Broad Institute)

**Our intention is to make GenBank intermediate files open and available for everyone to use, and to keep these files up-to-date.** The paths for specific files are the same under each endpoint, e.g. ``https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz``, ``s3://nextstrain-data/files/ncov/open/metadata.tsv.gz``, and ``gs://nextstrain-data/files/ncov/open/metadata.tsv.gz`` all exist. See below for a list of files that exist. If you're running workflows on AWS or GCP compute that fetch this data, please use the S3 or GS URLs, respectively, for cheaper (for us) and faster (for you) data transfers. Otherwise, please use the https://data.nextstrain.org URLs.

Note that even though the ``s3://nextstrain-data/`` and ``gs://nextstrain-data/`` buckets are public, the defaults for most S3 and GS clients require *some* user to be authenticated, though the specific user/account doesn't matter. In the rare case you need to access the S3 or GS buckets anonymously, the easiest way is to configure your inputs using ``https://nextstrain-data.s3.amazonaws.com/files/ncov/open/`` or ``https://storage.googleapis.com/nextstrain-data/files/ncov/open/`` URLs instead.

Depending on your execution environment, you may need to install additional Python dependencies for specific support of the different URL schemes (``https``, ``s3``, ``gs``). The workflow will produce an error at the start if additional dependencies are needed to fetch your configured inputs. Both ``https`` and ``s3`` should work out of the box in the standard Nextstrain Conda and Docker execution environments.

All available genomes and metadata
----------------------------------

Entire metadata & sequences data is uploaded from the ``ncov-ingest`` workflows for each of the ``gisaid`` and ``open`` sources:

-  ``metadata.tsv.gz`` and ``metadata.tsv.zst``
-  ``sequences.fasta.xz`` and ``sequences.fasta.zst``
-  ``nextclade.tsv.gz`` and ``nextclade.tsv.zst``
-  ``aligned.fasta.xz`` and ``aligned.fasta.zst`` (Alignment via `Nextclade <https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files.html#aligned-nucleotide-sequences>`__. The default reference genome is `MN908947 <https://www.ncbi.nlm.nih.gov/nuccore/MN908947>`__ (Wuhan-Hu-1))
-  ``additional_info.tsv.gz`` and ``additional_info.tsv.zst`` (GISAID only)
-  ``flagged_metadata.txt.gz`` and ``flagged_metadata.txt.zst`` (GISAID only)

Subsampled datasets
-------------------

Our GISAID and GenBank (open) profiles each define 7 builds (a Global build and one build per region: Africa, Asia, Europe, Oceania, North and South America). Each of these is a different subsample of the entire dataset, and each will result in the following intermediates uploaded:

-  ``{build_name}/sequences.fasta.xz``
-  ``{build_name}/metadata.tsv.xz``
-  ``{build_name}/aligned.fasta.xz``
-  ``{build_name}/{build_name}.json`` (the main Auspice dataset file)
-  ``{build_name}/{build_name}_tip-frequencies.json``
-  ``{build_name}/{build_name}_root-sequence.json``

--------------

.. _remote-inputs-open-files:

Summary of available GenBank (open) files
-----------------------------------------

Each regional build (``global``, ``africa``, ``asia``, ``europe``, ``north-america``, ``oceania`` and ``south-america``) contains a subsampled set of approximately 4000 sequences. They are a good starting point if you are seeking a representative sample of data. Where available, this table also provides the URL for the resulting Auspice visualisation of the data.

   Please note that these files are uploaded in two batches (see above for details). This means that the full GenBank metadata and sequences are typically updated a couple of hours before the more processed files.

.. warning::
  The zstandard (zstd) files are not yet supported as direct inputs for the pipeline.

+-----------------------+-----------------------+------------------------------------------------------------------------------+
| description           | type                  | address                                                                      |
+=======================+=======================+==============================================================================+
| Full GenBank data     | metadata              | https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz                  |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz               |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz                 |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | metadata (zstd)       | https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.zst          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences (zstd)      | https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.zst       |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned (zstd)        | https://data.nextstrain.org/files/ncov/open/global/aligned.fasta.zst         |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| Global sample         | metadata              | https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz           |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.xz        |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/global/aligned.fasta.xz          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/global                                      |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| Africa sample         | metadata              | https://data.nextstrain.org/files/ncov/open/africa/metadata.tsv.xz           |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/africa/sequences.fasta.xz        |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/africa/aligned.fasta.xz          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/africa                                      |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| Asia sample           | metadata              | https://data.nextstrain.org/files/ncov/open/asia/metadata.tsv.xz             |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/asia/sequences.fasta.xz          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/asia/aligned.fasta.xz            |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/asia                                        |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| Europe sample         | metadata              | https://data.nextstrain.org/files/ncov/open/europe/metadata.tsv.xz           |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/europe/sequences.fasta.xz        |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/europe/aligned.fasta.xz          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/europe                                      |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| North America sample  | metadata              | https://data.nextstrain.org/files/ncov/open/north-america/metadata.tsv.xz    |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/north-america/sequences.fasta.xz |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/north-america/aligned.fasta.xz   |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/north-america                               |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| Oceania sample        | metadata              | https://data.nextstrain.org/files/ncov/open/oceania/metadata.tsv.xz          |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/oceania/sequences.fasta.xz       |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/oceania/aligned.fasta.xz         |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/oceania                                     |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| South America sample  | metadata              | https://data.nextstrain.org/files/ncov/open/south-america/metadata.tsv.xz    |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/south-america/sequences.fasta.xz |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/south-america/aligned.fasta.xz   |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | auspice               | https://nextstrain.org/ncov/open/south-america                               |
+-----------------------+-----------------------+------------------------------------------------------------------------------+

--------------

Starting your build from these intermediates
--------------------------------------------

Each workflow defines one or more inputs in the workflow config file.

In the simplest form, an input specifies a local path to some metadata and sequences, like so:

.. code:: yaml

   inputs:
     - name: example-data
       metadata: data/example_metadata.tsv
       sequences: data/example_sequences.fasta

Using the above table, we can easily modify this to create a build which uses the global subsample of GenBank data:

.. code:: yaml

   inputs:
     - name: global-representative-genbank-sample
       metadata: https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.gz
       sequences: https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.gz

To avoid unnecessarily aligning these sequences, we can instead start from the aligned sequences, like so:

.. code:: yaml

   inputs:
     - name: global-representative-genbank-sample
       metadata: https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.gz
       aligned: https://data.nextstrain.org/files/ncov/open/global/aligned.fasta.gz

The following starting points are available:

-  replace ``sequences`` with ``aligned`` (skips alignment)

Compressed vs uncompressed starting points
------------------------------------------

The workflow supports compressed metadata and sequences for any input stage. Files may be compressed using ``xz`` (``.xz``) or ``gzip`` (``.gz``) compression.
