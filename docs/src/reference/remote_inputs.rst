Remote inputs
=============

This page provides an overview of intermediate files which Nextstrain produces via daily workflow runs. Where appropriate, these files can be starting points for the `ncov workflow <https://github.com/nextstrain/ncov/>`__ (discussed below).

We have two GitHub repositories which routinely upload files to `S3 buckets <https://aws.amazon.com/s3/>`__: `ncov-ingest <https://github.com/nextstrain/ncov-ingest/>`__ and `ncov <https://github.com/nextstrain/ncov/>`__. Each of those runs separate pipelines for GISAID and open (Genbank and Robert Koch Institute (RKI)) data sources; these pipelines start with data curation and QC steps and end with the phylogenetic analyses you can see on `nextstrain.org <https://nextstrain.org/sars-cov-2/>`__

The GISAID data is stored at ``s3://nextstrain-ncov-private`` and is not publicly available, in line with the GISAID Terms of Use (this is used internally by Nextstrain).

The open (GenBank and RKI) data is publicly available at three endpoints:

-  ``https://data.nextstrain.org/files/ncov/open/``
-  ``s3://nextstrain-data/files/ncov/open/``
-  ``gs://nextstrain-data/files/ncov/open/`` (mirrored daily from S3 by the Broad Institute)

**Our intention is to make GenBank and RKI intermediate files open and available for everyone to use, and to keep these files up-to-date.** The paths for specific files are the same under each endpoint, e.g. ``https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst``, ``s3://nextstrain-data/files/ncov/open/metadata.tsv.zst``, and ``gs://nextstrain-data/files/ncov/open/metadata.tsv.zst`` all exist. See below for a list of files that exist. If you're running workflows on AWS or GCP compute that fetch this data, please use the S3 or GS URLs, respectively, for cheaper (for us) and faster (for you) data transfers. Otherwise, please use the https://data.nextstrain.org URLs.

Note that even though the ``s3://nextstrain-data/`` and ``gs://nextstrain-data/`` buckets are public, the defaults for most S3 and GS clients require *some* user to be authenticated, though the specific user/account doesn't matter. In the rare case you need to access the S3 or GS buckets anonymously, the easiest way is to configure your inputs using ``https://nextstrain-data.s3.amazonaws.com/files/ncov/open/`` or ``https://storage.googleapis.com/nextstrain-data/files/ncov/open/`` URLs instead.

Depending on your execution environment, you may need to install additional Python dependencies for specific support of the different URL schemes (``https``, ``s3``, ``gs``). The workflow will produce an error at the start if additional dependencies are needed to fetch your configured inputs. Both ``https`` and ``s3`` should work out of the box in the standard Nextstrain Conda and Docker execution environments.

All available genomes and metadata
----------------------------------

Entire metadata & sequences data is uploaded from the ``ncov-ingest`` workflows for each of the ``gisaid`` and ``open`` sources:

-  ``metadata.tsv.zst`` and ``metadata.tsv.gz``
-  ``sequences.fasta.zst`` and ``sequences.fasta.xz``
-  ``nextclade.tsv.zst`` and ``nextclade.tsv.gz``
-  ``aligned.fasta.zst`` and ``aligned.fasta.xz`` (Alignment via `Nextclade <https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/02-nuc-alignment.html>`__. The default reference genome is `MN908947 <https://www.ncbi.nlm.nih.gov/nuccore/MN908947>`__ (Wuhan-Hu-1))
-  ``additional_info.tsv.zst`` and ``additional_info.tsv.gz`` (GISAID only)
-  ``flagged_metadata.txt.zst`` and ``flagged_metadata.txt.gz`` (GISAID only)

The files compressed with Zstandard (``.zst``) will generally be faster to download (i.e. smaller in size) and faster to decompress than those compressed with ``xz``.

Subsampled datasets
-------------------

Our GISAID and open profiles each define 7 builds (a Global build and one build per region: Africa, Asia, Europe, Oceania, North and South America). Each of these is a different subsample of the entire dataset, and each will result in the following intermediates uploaded:

-  ``{build_name}/sequences.fasta.xz``
-  ``{build_name}/metadata.tsv.xz``
-  ``{build_name}/aligned.fasta.xz``
-  ``{build_name}/{build_name}.json`` (the main Auspice dataset file)
-  ``{build_name}/{build_name}_tip-frequencies.json``
-  ``{build_name}/{build_name}_root-sequence.json``

100k Subsamples
---------------

We also produce a subsample of the entire open dataset of around 100,000 samples.
This is particularly useful for development purposes or to run builds locally as the file sizes are typically around 10Mb (metadata) and 20Mb (sequences).
The data is chosen by sampling 50,000 samples from the previous 12 months and 50,000 prior to that, and within each sample we group by year, month and country in an attempt at even sampling.

--------------

.. _remote-inputs-open-files:

Summary of available open files
-----------------------------------------

Each regional build (``global``, ``africa``, ``asia``, ``europe``, ``north-america``, ``oceania`` and ``south-america``) contains a subsampled set of approximately 4000 sequences. They are a good starting point if you are seeking a representative sample of data. Where available, this table also provides the URL for the resulting Auspice visualisation of the data.

   Please note that these files are uploaded in two batches (see above for details). This means that the full GenBank metadata and sequences are typically updated a couple of hours before the more processed files.

+-----------------------+-----------------------+------------------------------------------------------------------------------+
| description           | type                  | address                                                                      |
+=======================+=======================+==============================================================================+
| Full open data        | metadata              | https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst                 |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/sequences.fasta.zst              |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned               | https://data.nextstrain.org/files/ncov/open/aligned.fasta.zst                |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | metadata (gz)         | https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz                  |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences (xz)        | https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz               |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | aligned (xz)          | https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz                 |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
| 100k sample           | metadata              | https://data.nextstrain.org/files/ncov/open/100k/metadata.tsv.xz             |
+-----------------------+-----------------------+------------------------------------------------------------------------------+
|                       | sequences             | https://data.nextstrain.org/files/ncov/open/100k/sequences.fasta.xz          |
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
       metadata: https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz
       sequences: https://data.nextstrain.org/files/ncov/open/global/sequences.fasta.xz

To avoid unnecessarily aligning these sequences, we can instead start from the aligned sequences, like so:

.. code:: yaml

   inputs:
     - name: global-representative-genbank-sample
       metadata: https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz
       aligned: https://data.nextstrain.org/files/ncov/open/global/aligned.fasta.xz

The following starting points are available:

-  replace ``sequences`` with ``aligned`` (skips alignment)

Compressed vs uncompressed starting points
------------------------------------------

The workflow supports compressed metadata and sequences for any input stage. Files may be compressed using Zstandard (``.zst``), xz (``.xz``), or gzip (``.gz``) compression.

Data origins
------------------------------------------

The data for our open dataset comes from `NCBI Genbank <https://www.ncbi.nlm.nih.gov/>`__ (via API), and the `Robert Koch Institute (RKI) <https://github.com/robert-koch-institut/SARS-CoV-2-Sequenzdaten_aus_Deutschland>`__. Some UK metadata is augmented with data available from `COG-UK <https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/>`__ (via CLIMB).
