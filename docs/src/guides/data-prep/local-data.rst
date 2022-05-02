Prepare your own local data
===========================

Prepare your own local data for quality control prior to submission to a public database.

.. contents:: Table of Contents
   :local:

Formatting your sequence data
-----------------------------

The first 2 lines in ``data/sequences.fasta`` look like this:

::

   >Wuhan-Hu-1/2019
   ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATC.....

**The first line is the ``strain`` or ``name`` of the sequence.** Lines with names in FASTA files always start with the ``>`` character (this is not part of the name), and may not contain spaces or ``()[]{}|#><``. Note that “strain” here carries no biological or functional significance and should largely be thought of as synonymous with “sample.”

The sequence itself is a `consensus genome <https://en.wikipedia.org/wiki/Consensus_sequence#:~:text=In%20molecular%20biology%20and%20bioinformatics,position%20in%20a%20sequence%20alignment.>`__.

**By default, sequences less than 27,000 bases in length or with more than 3,000 ``N`` (unknown) bases are omitted from the analysis.** For a basic QC and preliminary analysis of your sequence data, you can use `clades.nextstrain.org <https://clades.nextstrain.org/>`__. This tool will check your sequences for excess divergence, clustered differences from the reference, and missing or ambiguous data. In addition, it will assign nextstrain clades and call mutations relative to the reference.

Formatting your metadata
------------------------

Nextstrain accommodates many kinds of metadata, so long as it is in a ``TSV`` format. A ``TSV`` is a text file, where each row (line) represents a sample and each column (separated by tabs) represents a field.

.. note::

   If you're unfamiliar with TSV files, don't fret; it's straightforward to export these directly from Excel, which we'll cover shortly.

Here's an example of the first few columns of the metadata for a single strain, including the header row. *(Spacing between columns here is adjusted for clarity, and only the first 6 columns are shown).*

::

   strain              virus  gisaid_epi_isl  genbank_accession   date        region   ...
   NewZealand/01/2020  ncov   EPI_ISL_413490  ?                   2020-02-27  Oceania  ...

:doc:`See the reference guide on metadata fields for more details <../../reference/metadata-fields>`.

Required metadata
~~~~~~~~~~~~~~~~~

A valid metadata file must include the following fields:

+------------------------+---------------------------------------------------------------------------------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------+
| Field                  | Example value                                                                         | Description                 | Formatting                                                                                                                    |
+========================+=======================================================================================+=============================+===============================================================================================================================+
| ``strain`` or ``name`` | ``Australia/NSW01/2020``                                                              | Sample name / ID            | Each header in the fasta file must exactly match a ``strain`` value in the metadata. Characters ``()[]{}|#><`` are disallowed |
+------------------------+---------------------------------------------------------------------------------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------+
| ``date``               | ``2020-02-27``, ``2020-02-XX``, ``2020-XX-XX``                                        | Date of *sampling*          | ``YYYY-MM-DD``; ambiguities can be indicated with ``XX``                                                                      |
+------------------------+---------------------------------------------------------------------------------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------+
| ``virus``              | ``ncov``                                                                              | Pathogen name               | Needs to be consistent                                                                                                        |
+------------------------+---------------------------------------------------------------------------------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------+
| ``region``             | ``Africa``, ``Asia``, ``Europe``, ``North America``, ``Oceania`` or ``South America`` | Global region of *sampling* |                                                                                                                               |
+------------------------+---------------------------------------------------------------------------------------+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------+

Please be aware that **our current workflow will filter out any genomes with an unknown date - you can change this in your own workflow.**

Missing metadata
~~~~~~~~~~~~~~~~

Missing data is to be expected for certain fields. In general, **missing data is represented by an empty string or a question mark character.** There is one important difference: if a discrete trait reconstruction (e.g. via ``augur traits``) is to be run on this column, then a value of ``?`` will be inferred, whereas the empty string will be treated as missing data in the output. See below for how to represent uncertainty in sample collection date.

General formatting tips
~~~~~~~~~~~~~~~~~~~~~~~

-  **The order of the fields doesn't matter**; but if you are going to join your metadata with the global collection then it's easiest to keep them in the same order!
-  **Not all fields are currently used**, but this may change in the future.
-  Data is **case sensitive**.
-  The **“geographic” columns, such as “region” and “country” will be used to plot the samples on the map**. Adding a new value to these columns isn't a problem at all, but there are a few extra steps to take; see :doc:`../workflow-config-file`.
-  **You can color by any of these fields in the Auspice visualization**. Which exact columns are used, and which colors are used for each value is completely customizable; see :doc:`../customizing-visualization`.

Formatting metadata in Excel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also create a TSV file in Excel. However, due to issues with auto-formatting of certain fields in Excel (like dates), we don't recommend this as a first option. If you do edit a file in Excel, open it afterwards in a text editor to check it looks as it should.

1. Create a spreadsheet where each row is a sample, and each column is a metadata field
2. Ensure your spreadsheet meets the requirements outlined above. Pay special attention to date formats; see `this guide to date formatting in Excel <https://support.microsoft.com/en-us/office/format-a-date-the-way-you-want-8e10019e-d5d8-47a1-ba95-db95123d273e?ui=en-us&rs=en-us&ad=us>`__.
3. Click on ``File > Save as``
4. Choose ``Text (Tab delimited) (*.txt)`` and enter a filename ending in ``.tsv``
