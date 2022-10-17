.. cssclass:: configuration-reference

Workflow config file reference
==============================

This is the detailed reference for sections in a :term:`workflow config file <config file>`. For example use cases, see :doc:`../guides/workflow-config-file`.

.. contents:: Table of Contents
   :local:
   :depth: 2

Primary configuration
^^^^^^^^^^^^^^^^^^^^^

Parameters in this section define the main inputs and outputs of the workflow, as well as the commonly used ``subsampling`` rule.
Often these will be the only parameters you need to modify.

inputs
------

-  type: array
-  description: A list of named input datasets to use for the workflow. Input order determines the precedence of genome sequences and metadata such that later datasets override earlier datasets. Each input must define a ``name``, a path to ``metadata``, and a path to sequences at one of many possible starting points. The workflow merged all input metadata and sequences into a single metadata and sequences file prior to subsampling.
-  required

   -  ``name``
   -  ``metadata``
   -  ``sequences`` or ``aligned``

-  examples:

.. code:: yaml

   inputs:
     - name: example-data
       metadata: data/example_metadata.tsv.xz
       sequences: data/example_sequences.fasta.xz
     - name: prealigned-data
       metadata: data/other_metadata.tsv.xz
       aligned: data/other_aligned.fasta.xz

Valid attributes for list entries in ``inputs``:

.. contents::
   :local:

name
~~~~

-  type: string
-  description: Name of the current dataset. Names cannot contain spaces, as they correspond to files on the file system.
-  examples:

   -  ``example-data``
   -  ``gisaid``
   -  ``washington``
   -  ``north-america``

metadata
~~~~~~~~

-  type: string
-  description: Path to a local or remote (S3, HTTP(S), GS) tab-delimited metadata file supported by Augur. Metadata can be uncompressed or compressed.
-  examples:

   -  ``data/example_metadata.tsv``
   -  ``data/example_metadata.tsv.xz``
   -  ``s3://your-bucket/metadata.tsv.gz``
   -  ``https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz``

sequences
~~~~~~~~~

-  type: string
-  description: Path to a local or remote (S3, HTTP(S), GS) FASTA file with \**_un_aligned*\* genome sequences. Sequences can be uncompressed or compressed.
-  examples:

   -  ``data/example_sequences.fasta``
   -  ``data/example_sequences.fasta.xz``
   -  ``s3://your-bucket/sequences.fasta.gz``
   -  ``https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz``

aligned
~~~~~~~

-  type: string
-  description: Path to a local or remote (S3, HTTP(S), GS) FASTA file with **aligned** genome sequences. Sequences can be uncompressed or compressed.
-  examples:

   -  ``data/aligned.fasta``
   -  ``data/aligned.fasta.xz``
   -  ``s3://your-bucket/aligned.fasta.gz``
   -  ``https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz``


builds
------

-  type: object
-  description: Named builds to produce by the workflow from the given inputs. Builds are indexed by name and include any number of build attributes that can be used to control subsampling, Auspice configuration, and more.

.. warning::

   Build names only allow alphanumeric characters, underscores, and hyphens (``A-Z``, ``a-z``, ``0-9``, ``_``, ``-``), but must not contain ``tip-frequencies`` as it is a special string used for Nextstrain builds.

   Note that these are not allowed:

   - Periods (``.``)

-  examples:

.. code:: yaml

   builds:

     # the following build (dataset) will include all samples provided in the inputs
     everything:
      subsampling_scheme: all

     # this will use a predefined subsampling scheme (see subsampling section for details)
     washington:
       region: North America
       country: USA
       division: Washington
       subsampling_scheme: country

     # this will use a custom subsampling scheme that you provide
     # which will have access to the provided `my_param`
     washington:
       subsampling_scheme: my_scheme
       my_param: some value

Valid attributes for entries in ``builds``:

.. contents::
   :local:

<user-specified key>
~~~~~~~~~~~~~~~~~~~~

-  description: Builds support any named attributes that can be referenced by subsampling schemes. In the above example, "country" and "division" are examples of this.


auspice_config
~~~~~~~~~~~~~~

-  type: string
-  description: Path to a build-specific Auspice configuration JSON.

colors
~~~~~~

-  type: string
-  description: Path to a build-specific color map to use in Auspice.

.. _configuration-builds-description:

description
~~~~~~~~~~~

-  type: string
-  description: Path to a build-specific Markdown file to use as a description in Auspice. Overrides :ref:`files.description <configuration-files-description>`.

region
~~~~~~

-  type: string
-  description: Name of the region the corresponding build belongs to (based on standard values in the ``region`` metadata field).

.. warning::

   The presence of a ``region`` key will result in the metadata being adjusted in potentially surprising ways.
   For all metadata rows that are not in this region, ``location`` will be removed (set to an empty string), and ``division`` and ``country`` will be changed to their corresponding region.
   Additionally, a ``focal`` column will be added, with True/False values depending on if the row matches the provided region.

subclades
~~~~~~~~~

-  type: string
-  description: Path to a build-specific `Augur clade definition file <https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations>`__ to combine with the curated clades defined by ``files: clades``.

subsampling_scheme
~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Name of the subsampling scheme defined in ``subsampling`` to use for the current build.
-  default: ``"all"``. In practice, this means that no subsampling will be performed.

title
~~~~~

-  type: string
-  description: Build-specific title to provide to ``augur export`` and display as the title of the analysis in Auspice.


.. _configuration-subsampling:

subsampling
-----------

-  type: object
-  description: Schemes for subsampling data prior to phylogenetic inference to avoid sampling bias or focus an analysis on specific spatial and/or temporal scales. `See the SARS-CoV-2 tutorial for more details on defining subsampling schemes <../reference/customizing-analysis.html#subsampling>`__.

Predefined subsampling schemes are:

- ``all``
- ``region``
- ``region_global``
- ``region_grouped_by_country``
- ``country``
- ``division``
- ``location``

See `defaults/parameters.yaml <https://github.com/nextstrain/ncov/blob/master/defaults/parameters.yaml>`__ for definitions.

Each named subsampling scheme supports the following attributes that the workflow passes to ``augur filter``.

.. contents::
   :local:

group_by
~~~~~~~~

-  type: string
-  description: Space-delimited list of metadata columns to group records by prior to subsampling to the requested or calculated number of sequences per group.
-  examples:

   -  ``year month``
   -  ``region year month``

seq_per_group
~~~~~~~~~~~~~

-  type: integer
-  description: Number of sequences to select per group of records in groups specified by ``group_by``. The total number of sequences selected for each subsampling rule will be no more than the number of groups times this number of sequences per group. This parameter must be used with the ``group_by`` parameter.

max_sequences
~~~~~~~~~~~~~

-  type: integer
-  description: Maximum number of sequences to select for the current subsampling rule. When used with the ``group_by`` parameter, Augur will calculate the number of sequences per group. When used without the ``group_by`` parameter, Augur will select this number of sequences at random from all available sequences. When probabilistic sampling is enabled by the ``sampling_scheme`` parameter, the total number of strains actually selected will be more or less than this value due to the underlying Poisson sampling process.

sampling_scheme
~~~~~~~~~~~~~~~

-  type: string
-  description: A flag to pass to ``augur filter`` that specifies whether to enable probabilistic sampling or not. Probabilistic sampling is useful when there are more groups than requested sequences.
-  default: ``--probabilistic-sampling`` (Augur's default)
-  examples:

   -  ``--probabilistic-sampling``
   -  ``--no-probabilistic-sampling``

.. _exclude-1:

exclude
~~~~~~~

-  type: string
-  description: Argument to pass to ``augur filter`` to exclude records based on specific values in metadata columns. This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below.
-  examples:

   -  ``"--exclude-where 'region!=Africa'"``
   -  ``"--exclude-where 'region!={region}'"``

.. _include-1:

include
~~~~~~~

-  type: string
-  description: Argument to pass to ``augur filter`` to include records based on specific values in metadata columns regardless of other filters applied during subsampling (i.e., strains for which the include test evaluates to true will always be included if they exist in the metadata and sequences). This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below.
-  examples:

   -  ``--include-where 'region=Africa'``
   -  ``--include-where 'region={region}'``

query
~~~~~

-  type: string
-  description: Argument to pass to ``augur filter`` to select specific records by testing values in metadata columns. This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below. Query values support `pandas Dataframe query syntax <https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query>`__ treating the metadata as a data frame.
-  examples:

   -  ``--query "division == 'Washington'"``
   -  ``--query "division == '{division}'"``
   -  ``--query "(country == '{country}') & (division == '{division}')"``
   -  ``--query "division != '{division}'"``

.. _exclude_ambiguous_dates_by-1:

exclude_ambiguous_dates_by
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Level date ambiguity used to exclude strains from the analysis by ``augur filter --exclude-ambiguous-dates-by``
-  examples:

   -  ``any``
   -  ``day``
   -  ``month``
   -  ``year``

.. _min_date-2:

min_date
~~~~~~~~

-  type: string
-  description: Argument to ``augur filter`` to set the minimum collection date for strains to include in the subsampling set. See :doc:`augur filter docs <augur:usage/cli/filter>` for supported date formats.
-  examples:

   -  ``--min-date 2019-10-01``
   -  ``--min-date 2019.74``

.. _max_date-1:

max_date
~~~~~~~~

-  type: string
-  description: Argument to ``augur filter`` to set the maximum collection date for strains to include in the subsampling set. See :doc:`augur filter docs <augur:usage/cli/filter>` for supported date formats.
-  examples:

   -  ``--max-date 2021-04-01``
   -  ``--max-date 2021.25``

priorities
~~~~~~~~~~

-  type: object
-  description: Parameters to prioritize strains selected for the current subsampling rule. Currently, the workflow supports two ``type``\ s of priority, ``proximity`` and ``file``.
-  description [proximity]: ``proximity`` selects samples that are genetically similar to the ``focus`` sample set; the ``focus`` sample set must be a rule in the current subsampling scheme.
-  example [proximity]:

.. code:: yaml

   subsampling:
     my-scheme:
       my-first-rule:
         max_sequences: 10
       my-second-rule:
         max_sequences: 10
         # Prioritize sequences that are genetically similar to
         # sequences in the sequences selected by the
         # `my-first-rule` rule.
         priorities:
           type: proximity
           focus: my-first-rule

-  description [file]: ``file`` selects samples based on arbitrarily-defined rankings in a TSV file formatted as ``strain\tnumber``. The numbers are only used to sort the samples, and are therefore arbitrary. Higher values = higher priority.

-  example [file]:

.. code:: yaml

   subsampling:
     my-scheme:
       my-first-rule:
         max_sequences: 10
         group_by: "country"
         priorities:
           type: "file"
           file: "path/to/priorities.tsv"

::

   hCoV-19/USA/CZB-1234/2021   8.2
   hCoV-19/USA/CZB-2345/2021   0
   hCoV-19/USA/CZB-3456/2021   -3.1


Secondary configuration
^^^^^^^^^^^^^^^^^^^^^^^

These parameters are other high-level parameters which may affect multiple Snakemake rules, or modify which rules are run.

default_build_name
------------------

-  type: string
-  description: Name to assign the default build when a user has not defined any other entries in the ``builds`` config.
-  default: ``default-build``


strip_strain_prefixes
---------------------

-  type: array
-  description: A list of prefixes to strip from strain names in metadata and sequence records to maintain consistent strain names when analyzing data from multiple sources.
-  default: ``["hCoV-19/", "SARS-CoV-2/"]``


auspice_json_prefix
-------------------

-  type: string
-  description: Prefix to use for Auspice JSON outputs. Change this value to produce JSONs named like ``auspice/<your_prefix>_global.json`` for a build named ``global``, for example. If you are using :doc:`Nextstrain's Community Sharing <docs.nextstrain.org:guides/share/community-builds>` to view your builds, set this value to your GitHub repository name and the ``ncov`` default. For example, if your repository is named ``evolution``, set ``auspice_json_prefix: evolution_ncov`` to get JSONs you can view your ``global`` build at https://nextstrain.org/community/*your_github_organization*/evolution/ncov/global.
-  default: ``ncov``


include_hcov19_prefix
---------------------

-  type: boolean
-  description: Prepend strain names with ``hCoV-19/`` per GISAID requirements for web display
-  default: ``false``


title
-----

-  type: string
-  description: Title to provide to ``augur export`` and display as the title of the analysis in Auspice. Note that this is only used if a title is not defined for the individual build in the ``builds`` object.


genes
-----

-  type: array
-  description: A list of genes for which ``nextalign`` should generate amino acid sequences during the alignment process. Gene names must match the names provided in the gene map from the ``annotation`` parameter.
-  default: ``["ORF1a", "ORF1b", "S", "ORF3a", "M", "N"]``
-  used in rules: ``align``, ``build_align``, ``translate``, ``mutational_fitness``


active_builds
-------------

-  type: string
-  description: Comma-delimited list of names of builds to run (allowing a subset of all builds to be specified). You only need to use this parameter if you want to run a subset of the builds defined in ``builds``.
-  examples

   -  ``global``
   -  ``global,africa,north-america``



conda_environment
-----------------

-  type: string
-  description: Path to a Conda environment file to use for the workflow when the workflow is run with `Snakemake's ``--use-conda`` flag <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management>`__.
-  default: ``workflow/envs/nextstrain.yaml``

custom_rules
------------

-  type: array
-  description: List of paths to Snakemake files to include in the workflow, allowing users to inject their own rules at the beginning or the end of the workflow (e.g., to pre-process data prior to the workflow, annotate outputs from the workflow, etc.).
-  examples

   -  ``- workflow/snakemake_rules/export_for_nextstrain.smk``
   -  ``- nextstrain_profiles/nextstrain-gisaid/subsampling_ranges.smk``


localrules
----------

-  type: string
-  description: Path to a Snakemake file to include in the workflow. This parameter is redundant with ``custom_rules`` and may be deprecated soon.




files
-----

-  type: object
-  description: Additional files used to configure tools used by the workflow (e.g., alignment references, names of strains to exclude during filtering, etc.).
- Valid attributes:

.. contents::
   :local:

include
~~~~~~~

-  type: string
-  description: Path to a file with list of strains (one name per line) to include in the analysis regardless of priorities or subsampling during filtering.
-  default: ``defaults/include.txt``
-  used in rules: ``subsample``, ``filter``

exclude
~~~~~~~

-  type: string
-  description: Path to a file with list of strains (one name per line) to exclude from the analysis.
-  default: ``defaults/exclude.txt``
-  used in rules: ``subsample``, ``filter``

reference
~~~~~~~~~

-  type: string
-  description: Path to a GenBank-formatted sequence to use for sequence translation
-  default: ``defaults/reference_seq.gb``
-  used in rules: ``translate``

alignment_reference
~~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Path to a FASTA-formatted sequence to use for alignment with ``nextalign``
-  default: ``defaults/reference_seq.fasta``
-  used in rules: ``align``, ``proximity_score`` (subsampling), ``build_align``, ``build_mutation_summary``

annotation
~~~~~~~~~~

-  type: string
-  description: Path to a GFF-formatted annotation of gene coordinates (e.g., a “gene map”) for use by ``nextalign`` and mutation summaries.
-  default: ``defaults/annotation.gff``
-  used in rules: ``align``, ``build_align``, ``build_mutation_summary``

outgroup
~~~~~~~~

-  type: string
-  description: No longer used.

ordering
~~~~~~~~

-  type: string
-  description: Path to tab-delimited mapping of metadata attributes (first column) to corresponding values (second column) with rows ordered by the desired appearance in the Nextstrain color legend. This mapping and ordering is manually curated by the Nextstrain team and updates regularly. Along with the ``color_schemes`` file, this file is used to generate a build-specific color map for use by Auspice.
-  default: ``defaults/color_ordering.tsv``
-  used in rules: ``colors``

color_schemes
~~~~~~~~~~~~~

-  type: string
-  description: Path to a list of tab-delimited and manually curated categorical color schemes for N total categories where row one defines one color, row two define two colors, and so on. Along with the ``ordering`` file, this file is used to generate a build-specific color map for use by Auspice.
-  default: ``defaults/color_schemes.tsv``
-  used in rules: ``colors``

.. _auspice_config-1:

auspice_config
~~~~~~~~~~~~~~

-  type: string
-  description: Path to an Auspice configuration JSON file used by ``augur export``. Note that this is only used if a build does not define its own ``auspice_config`` (in the ``builds`` config section).
-  default: ``defaults/auspice_config.json``
-  used in rules: ``export``

lat_longs
~~~~~~~~~

-  type: string
-  description: Path to a tab-delimited mapping of geographic scales (e.g., ``location`` ,\ ``division``, etc.), geographic names (e.g., ``King County``), and corresponding latitude and longitude values for the given place name. This mapping is manually curated by the Nextstrain team and updates regularly.
-  default: ``defaults/lat_longs.tsv``
-  used in rules: ``export``

.. _configuration-files-description:

description
~~~~~~~~~~~

-  type: string
-  description: Path to a Markdown file to use as a description in Auspice for all builds. Overridden per-build by :ref:`builds.description <configuration-builds-description>`.
-  default: ``defaults/description.md``
-  used in rules: ``export``

clades
~~~~~~

-  type: string
-  description: Path to `an Augur clade definition file <https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations>`__ where each row is a tab-delimited mapping of clade name to a gene, site (i.e., position), and alternate allele at that site for the corresponding clade.
-  default: ``defaults/clades.tsv``
-  used in rules: ``emerging_lineages``, ``clades``

emerging_lineages
~~~~~~~~~~~~~~~~~

-  type: string
-  description: Path to `an Augur clade definition file <https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations>`__ for emerging lineages of concern that may be a subset or variation of the lineages defined by the ``clades`` parameter or Pangolin lineages.
-  default: ``defaults/emerging_lineages.tsv``
-  used in rules: ``emerging_lineages``


Per-Rule configuration
^^^^^^^^^^^^^^^^^^^^^^

Each top-level parameter here corresponds to a single Snakemake rule.
Note that ``subsampling`` is a commonly used rule configuration which is described separately in the Primary configuration section.

sanitize_metadata
-----------------

-  type: object
-  description: Parameters to configure how to sanitize metadata to a Nextstrain-compatible format. The sanitize metadata script resolves duplicate records using database ids, parses a GISAID-style location field into Nextstrain-style location fields, strips prefixes from strain names, and renames fields in that order.
- Valid attributes:

.. contents::
   :local:

metadata_id_columns
~~~~~~~~~~~~~~~~~~~

-  type: object
-  description: A list of valid strain name columns in the metadata. The sanitize metadata script will check attempt to use the first of these columns that exists in the metadata. It will exit with an error, if none of the columns exist.
-  default:

.. code:: yaml

     - strain
     - name
     - "Virus name"

database_id_columns
~~~~~~~~~~~~~~~~~~~

-  type: object
-  description: A list of columns representing external database ids for metadata records. These unique ids represent a snapshot of data at a specific time for a given strain name. The sanitize metadata script resolves duplicate metadata records for the same strain name by selecting the record with the latest database id. Multiple database id columns allow the script to resolve duplicates when one or more columns has ambiguous values (e.g., “?”). Deduplication occurs before renaming of columns, so the default values include GISAID's own “Accession ID” as well as Nextstrain-style database ids.
-  default:

.. code:: yaml

     - "Accession ID"
     - gisaid_epi_isl
     - genbank_accession

error_on_duplicate_strains
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  type: boolean
-  description: Exit the sanitize metadata script with an error when any strains have multiple records in the metadata. The script writes list of all duplicate strains to a file named like ``<input>.duplicates.txt`` that users can review and use to address unexpected duplicates.
-  default: ``false``

parse_location_field
~~~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Field in the metadata that stores GISAID-formatted location details (e.g., ``North America / USA / Washington``) to be parsed into ``region``, ``country``, ``division``, and ``location`` fields.
-  default: ``Location``

rename_fields
~~~~~~~~~~~~~

-  type: array
-  description: List of key/value pairs mapping fields in the input metadata to rename to another value in the sanitized metadata.
-  default:

.. code:: yaml

       - "Virus name=strain"
       - "Type=type"
       - "Accession ID=gisaid_epi_isl"
       - "Collection date=date"
       - "Additional location information=additional_location_information"
       - "Sequence length=length"
       - "Host=host"
       - "Patient age=patient_age"
       - "Gender=sex"
       - "Clade=GISAID_clade"
       - "Pango lineage=pango_lineage"
       - "Pangolin version=pangolin_version"
       - "Variant=variant"
       - "AA Substitutions=aa_substitutions"
       - "aaSubtitutions=aa_substitutions"
       - "Submission date=date_submitted"
       - "Is reference?=is_reference"
       - "Is complete?=is_complete"
       - "Is high coverage?=is_high_coverage"
       - "Is low coverage?=is_low_coverage"
       - "N-Content=n_content"
       - "GC-Content=gc_content"



combine_sequences_for_subsampling
---------------------------------

-  type: object
-  description: Configuration of logic to combine sequences from multiple input files into a single file for subsampling.
- Valid attributes:

.. contents::
   :local:


warn_about_duplicates
~~~~~~~~~~~~~~~~~~~~~

-  type: boolean
-  description: Warn users about duplicate sequences identified when merging input sequences and print a list of duplicates to standard out (and log files). Set this to ``false`` to get an error and stop the workflow when duplicates are detected.
-  default: ``true``


priorities
----------

-  type: object
-  description: Configures how proximities are calculated, which is used by subsampling schemes which specify it.
- Valid attributes:

.. contents::
   :local:

crowding_penalty
~~~~~~~~~~~~~~~~

-  type: float
-  description: used when calculating ``priority scores`` during subsampling to decrease the number of identical samples that are included in the tree during random subsampling to provide a broader picture of the viral diversity in your dataset.
-  examples:

.. code:: yaml

   priorities:
     crowding_penalty: 0.0
     # You may wish to set `crowding_penalty = 0.0` (default value = `0.1`) if you are interested in seeing as many samples as possible that are closely related to your `focal` set.

.. _title-1:

run_pangolin
------------

-  type: boolean
-  description: Enable annotation of Pangolin lineages for a given build's subsampled sequences.
-  default: ``false``

.. _workflow-config-mask:

mask
----

-  type: object
-  description: Parameters for masking of invalid or problematic nucleotides in aligned sequences. In addition to the configurable parameters below, the workflow also always masks terminal gaps in the given alignment.
-  Valid attributes:

.. contents::
   :local:

mask_from_beginning
~~~~~~~~~~~~~~~~~~~

-  type: integer
-  description: Number of bases to mask from the beginning alignment.
-  default: ``100``

mask_from_end
~~~~~~~~~~~~~

-  type: integer
-  description: Number of bases to mask from the end alignment.
-  default: ``50``

mask_sites
~~~~~~~~~~

-  type: string
-  description: Space-delimited string of 1-based genomic sites to mask
-  default: ``"13402 24389 24390"``



.. _workflow-config-filter:

filter
------

-  type: object
-  description: Filters to apply to strain metadata and sequences prior to subsampling and tree inference. The workflow applies an implicit filter on the maximum collection dates later than today.
- Valid attributes:

.. contents::
   :local:

min_length
~~~~~~~~~~

-  type: integer
-  description: Minimum number of valid nucleotides (A, C, T, or G) for a genome to be included in the analysis by ``augur filter --min-length``.
-  default: ``27000``

.. note::
   The ``min_length`` filter is applied to the **masked** sequences, not the original input sequences.
   Depending on your :ref:`mask config parameters <workflow-config-mask>`, the masked sequences may contain more Ns than the original sequences and fail the ``min_length`` filter.

exclude_where
~~~~~~~~~~~~~

-  type: string
-  description: Conditional tests of metadata columns used to exclude strains from the analysis by ``augur filter --exclude-where``
-  default: ``"division='USA'"``

exclude_ambiguous_dates_by
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Level date ambiguity used to exclude strains from the analysis by ``augur filter --exclude-ambiguous-dates-by``
-  default: ``any``
-  examples:

   -  ``any``
   -  ``day``
   -  ``month``
   -  ``year``

min_date
~~~~~~~~

-  type: float or string
-  description: Minimum collection date for strains to include in the analysis used by ``augur filter --min-date``. See :doc:`augur filter docs <augur:usage/cli/filter>` for supported date formats.
-  default: ``2019.74``

skip_diagnostics
~~~~~~~~~~~~~~~~

-  type: boolean
-  description: Skip filtering by Nextclade quality control metrics like clock rate deviation, number of SNP clusters, possible contaminations, etc.
-  default: ``false``



tree
----

-  type: object
-  description: Parameters for phylogenetic inference by ``augur tree``. The tree “method” is hardcoded to ``iqtree``.
-  Valid attributes:

.. contents::
   :local:

tree-builder-args
~~~~~~~~~~~~~~~~~

-  type: string
-  description: Arguments specific to the tree method (``iqtree``) to be passed through to the tree builder command run by ``augur tree``.
-  default: ``'-ninit 10 -n 4'``



refine
------

-  type: object
-  description: Parameters for inference of time trees with ``augur refine``.
-  Valid attributes:

.. contents::
   :local:

root
~~~~

-  type: string
-  description: Rooting mechanism or strain name(s) whose sequences should be used to root the time tree. Only one or two (space-delimited) strain names are supported.
-  default: ``Wuhan/Hu-1/2019``
-  examples:

   -  ``best``
   -  ``least-squares``
   -  ``min_dev``
   -  ``oldest``
   -  ``Wuhan/Hu-1/2019``

clock_rate
~~~~~~~~~~

-  type: float
-  description: Fixed clock rate to use for time tree calculations.
-  default: ``0.0008``

clock_std_dev
~~~~~~~~~~~~~

-  type: float
-  description: Standard deviation of the fixed ``clock_rate`` estimate.
-  default: ``0.0004``

coalescent
~~~~~~~~~~

-  type: float or string
-  description: Coalescent timescale in units of inverse clock rate (float), optimized as a scalar (“opt”), or skyline (“skyline”).
-  default: ``skyline``
-  examples:

   -  ``opt``
   -  ``skyline``

date_inference
~~~~~~~~~~~~~~

-  type: string
-  description: Assign internal nodes to their jointly or marginally most likely dates.
-  default: ``marginal``
-  examples:

   -  ``marginal``
   -  ``joint``

divergence_unit
~~~~~~~~~~~~~~~

-  type: string
-  description: Units in which sequence divergence is reported.
-  default: ``mutations``
-  examples:

   -  ``mutations``
   -  ``mutations-per-site``

clock_filter_iqd
~~~~~~~~~~~~~~~~

-  type: integer
-  description: Remove tips that deviate more than this number of interquartile ranges from the root-to-tip by time regression. Disable clock filtering by specifying ``0``
-  default: ``8``

keep_polytomies
~~~~~~~~~~~~~~~

-  type: boolean
-  description: Do not attempt to resolve polytomies.
-  default: ``false``

no_timetree
~~~~~~~~~~~

-  type: boolean
-  description: Do not produce a time tree.
-  default: ``false``


traits
------

-  type: object
-  description: Parameters for inference of ancestral traits by ``augur traits`` with support for default traits and build-specific traits.
-  examples:

.. code:: yaml

   traits:
     default:
       sampling_bias_correction: 2.5
       columns: ["country"]
     washington:
       # Override default sampling bias correction for
       # "washington" build and continue to use default
       # trait columns.
       sampling_bias_correction: 5.0

Each named traits configuration (``default`` or build-named) supports the following attributes:

.. contents::
   :local:

sampling_bias_correction
~~~~~~~~~~~~~~~~~~~~~~~~

-  type: float
-  description: A rough estimate of how many more events would have been observed if sequences represented an even sample. :doc:`See the documentation for augur traits for more details <augur:usage/cli/traits>`.
-  default: ``2.5``

columns
~~~~~~~

-  type: array
-  description: A list of columns from the metadata for which ancestral trait values should be inferred for ancestral nodes.
-  default: ``["country"]``


frequencies
-----------
- Valid attributes:

.. contents::
   :local:

.. _min_date-1:

min_date
~~~~~~~~

-  type: float or string
-  description: Earliest date to estimate frequencies for. See :doc:`augur filter docs <augur:usage/cli/filter>` for supported date formats.
-  default: without value supplied, defaults to 1 year before present

max_date
~~~~~~~~

-  type: float or string
-  description: Earliest date to estimate frequencies for. Specifying ``max_date`` overrides ``recent_days_to_censor``. See :doc:`augur filter docs <augur:usage/cli/filter>` for supported date formats.
-  default: without value supplied, defaults to today's date minus ``recent_days_to_censor`` parameter

recent_days_to_censor
~~~~~~~~~~~~~~~~~~~~~

-  type: integer
-  description: How many days back from today's date should samples be hidden from frequencies calculations? This is in place to help with sampling bias where some regions have faster sequencing turnarounds than other regions.
-  default: without value supplied, defaults to ``0``

pivot_interval
~~~~~~~~~~~~~~

-  type: integer
-  description: Number of units between frequency estimates based on the units defined in the ``pivot_interval_units`` parameter. A “pivot” corresponds to a time point when frequencies are estimated.
-  default: ``1``

pivot_interval_units
~~~~~~~~~~~~~~~~~~~~

-  type: string
-  description: Unit of pivot interval spacing for frequency estimation.
-  default: ``weeks``
-  examples:

   -  ``weeks``
   -  ``months``

narrow_bandwidth
~~~~~~~~~~~~~~~~

-  type: float
-  description: Variance of the KDE normal distribution in numeric floating point years (e.g., one month ~= 30 days ~= 0.08 years). This bandwidth value controls the smoothing of frequency estimates with higher values producing smoother estimates.
-  default: ``0.05``

proportion_wide
~~~~~~~~~~~~~~~

-  type: float
-  description: Proportion of a second KDE normal distribution to add to each initial normal distribution already parameterized by the ``narrow_bandwidth`` parameter.
-  default: ``0.0``

minimal_frequency
~~~~~~~~~~~~~~~~~

-  Unused

stiffness
~~~~~~~~~

-  Unused

inertia
~~~~~~~

-  Unused


logistic_growth
---------------

-  type: object
-  description: Parameters for estimation of logistic clade growth based on logit-transformed clade frequencies.
-  Valid attributes:

.. contents::
   :local:

delta_pivots
~~~~~~~~~~~~

-  type: integer
-  description: Calculate logistic growth over the last N pivots which corresponds to N times the amount of time represented by the ``pivot_interval_units`` in the frequencies configuration.
-  default: ``6``

min_tips
~~~~~~~~

-  type: integer
-  description: The minimum number of tips a clade must have before its logistic growth is calculated.
-  default: ``50``

min_frequency
~~~~~~~~~~~~~

-  type: float
-  description: The minimum current frequency for a clade to have its logistic growth calculated.
-  default: ``0.000001``

max_frequency
~~~~~~~~~~~~~

-  type: float
-  description: The maximum current frequency for a clade to have its logistic growth calculated.
-  default: ``0.95``



ancestral
---------

-  type: object
-  description: Configuration of augur ancestral command that infers ancestral sequences based on a tree.
- Valid attributes:

.. contents::
   :local:

inference
~~~~~~~~~

-  type: string
-  description: Calculate joint or marginal maximum likelihood ancestral sequence states
-  examples

   -  ``joint``
   -  ``marginal``


cluster
-------

-  type: object
-  description: Parameters for clustering of closely related strains
- Valid attributes:

.. contents::
   :local:

.. _min_tips-1:

min_tips
~~~~~~~~

-  type: integer
-  description: Number of tips to require in a polytomy to be considered part of a cluster.
-  default: ``3``

.. _group_by-1:

group_by
~~~~~~~~

-  type: string
-  description: Metadata column whose values should be used to determine whether closely related strains should be assigned to the same cluster. For example, the default column ensures that strains belong to the same division to be considered part of the same cluster.
-  default: ``division``

Internal, Nextstrain-only configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You shouldn't need to use parameters in this section unless you are running the core-Nextstrain builds. In addition to these parameters you will need AWS credentials set, e.g. in environment variables.

S3_DST_BUCKET
-------------

-  type: string
-  description: S3 bucket to store files from the ``upload`` rule in ``export_for_nextstrain.smk``. Currently only available to Nextstrain builds.

S3_DST_COMPRESSION
------------------

-  type: string
-  description: Compression format to use for files uploaded to S3 by the ``upload`` rule.
-  examples

   -  ``xz``
   -  ``gz``

S3_DST_ORIGINS
--------------

-  type: array
-  items:

   -  type: string

-  description: List of input names (i.e., “origins”) for which intermediate files should be uploaded to S3 by the ``upload`` rule.
-  examples

   -  ``["gisaid"]``

deploy_url
----------

-  type: string
-  description: URL to an S3 bucket where Auspice JSONs should be uploaded by the ``deploy`` rule of the Nextstrain workflows. Only valid for Nextstrain builds.

slack_channel
-------------

-  type: string
-  description: Slack channel to notify when Nextstrain builds start, fail, or get deployed. Only valid for Nextstrain builds.

slack_token
-----------

-  type: string
-  description: `Slack authentication token <https://api.slack.com/authentication/token-types>`__ required for the Slack API calls to notify the defined ``slack_channel``. Only valid for Nextstrain builds.

Unused Parameters
^^^^^^^^^^^^^^^^^

Documented here for completeness / historical accuracy.

partition_sequences
-------------------

-  Unused

reference_node_name
-------------------

-  Unused
