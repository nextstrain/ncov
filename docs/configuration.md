# Configuration parameters for Nextstrain SARS-CoV-2 workflow

## S3_DST_BUCKET
* type: string
* description: S3 bucket to store files from the `upload` rule in `export_for_nextstrain.smk`. Currently only available to Nextstrain builds.

## S3_DST_COMPRESSION
* type: string
* description: Compression format to use for files uploaded to S3 by the `upload`  rule.
* examples
	* `xz`
	* `gz`

## S3_DST_ORIGINS
* type: array
* items:
	* type: string
* description: List of input names (i.e., “origins”) for which intermediate files should be uploaded to S3 by the `upload` rule.
* examples
	* `["gisaid"]`

## active_builds
* type: string
* description: Comma-delimited list of names of builds to run (allowing a subset of all builds to be specified).
* examples
	* `global`
	* `global,africa,north-america`

## ancestral
* type: object
* description: Configuration of augur ancestral command that infers ancestral sequences based on a tree.

### inference
* type: string
* description: Calculate joint or marginal maximum likelihood ancestral sequence states
* examples
	* `joint`
	* `marginal`

## builds
* type: object
* description: Named builds to produce by the workflow from the given inputs. Builds are indexed by name and include any number of build attributes that can be used to control subsampling, Auspice configuration, and more.
* examples:
```yaml
builds:
  global:
    region: global
    subsampling_scheme: global

  washington:
    region: North America
    country: USA
    division: Washington
    subsampling_scheme: all
```
* required:
	* `region` (required to adjust regional metadata)

Builds support any named attributes that can be referenced by subsampling schemes. Builds also support the following specific attributes.

### auspice_config
* type: string
* description: Path to a build-specific Auspice configuration JSON.

### colors
* type: string
* description: Path to a build-specific color map to use in Auspice.

### description
* type: string
* description: Path to a build-specific Markdown file to use as a description in Auspice.

### region
* type: string
* description: Name of the region the corresponding build belongs to (based on standard values in the `region` metadata field).

### subclades
* type: string
* description: Path to a build-specific [Augur clade definition file](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations) to combine with the curated clades defined by `files: clades`.

### subsampling_scheme
* type: string
* description: Name of the subsampling scheme defined in `subsampling` to use for the current build.

### title
* type: string
* description: Build-specific title to provide to `augur export` and display as the title of the analysis in Auspice.

## combine_sequences_for_subsampling
* type: object
* description: Configuration of logic to combine sequences from multiple input files into a single file for subsampling.

### warn_about_duplicates
* type: boolean
* description: Warn users about duplicate sequences identified when merging input sequences and print a list of duplicates to standard out (and log files). Set this to `false` to get an error and stop the workflow when duplicates are detected.
* default: `true`

## conda_environment
* type: string
* description: Path to a Conda environment file to use for the workflow when the workflow is run with [Snakemake’s `--use-conda` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management).
* default: `workflow/envs/nextstrain.yaml`

## custom_rules
* type: array
* description:  List of paths to Snakemake files to include in the workflow, allowing users to inject their own rules at the beginning or the end of the workflow (e.g., to pre-process data prior to the workflow, annotate outputs from the workflow, etc.).
* examples
	* `- workflow/snakemake_rules/export_for_nextstrain.smk`
	* `- nextstrain_profiles/nextstrain-gisaid/subsampling_ranges.smk`

## default_build_name
* type: string
* description: Name to assign the default build when a user has not defined any other entries in the `builds` config.
* default: `default-build`

## exposure
* type: object
* description: Build-specific exposure history inference.

### default
* type: object
* description: Default exposure history inference at the country level.

## files
* type: object
* description: Additional files used to configure tools used by the workflow (e.g., alignment references, names of strains to exclude during filtering, etc.).

### include
* type: string
* description: Path to a file with list of strains (one name per line) to include in the analysis regardless of priorities or subsampling during filtering.
* default: `defaults/include.txt`

### exclude
* type: string
* description: Path to a file with list of strains (one name per line) to exclude from the analysis.
* default: `defaults/exclude.txt`

### reference
* type: string
* description: Path to a GenBank-formatted sequence to use for multiple sequence alignment with `augur align`
* default: `defaults/reference_seq.gb`

### alignment_reference
* type: string
* description: Path to a FASTA-formatted sequence to use for alignment with `nextalign` or `mafft`’s reference-based alignment
* default: `defaults/reference_seq.fasta`

### annotation
* type: string
* description: Path to a GFF-formated annotation of gene coordinates (e.g., a “gene map”) for use by `nextalign` and mutation summaries.
* default: `defaults/annotation.gff`

### outgroup
* type: string
* description: No longer used.
* default: `defaults/outgroup.fasta`

### ordering
* type: string
* description: Path to tab-delimited mapping of metadata attributes (first column) to corresponding values (second column) with rows ordered by the desired appearance in the Nextstrain color legend. This mapping and ordering is manually curated by the Nextstrain team and updates regularly. Along with the `color_schemes` file, this file is used to generate a build-specific color map for use by Auspice.
* default: `defaults/color_ordering.tsv`

### color_schemes
* type: string
* description: Path to a list of tab-delimited and manually curated categorical color schemes for N total categories where row one defines one color, row two define two colors, and so on. Along with the `ordering` file, this file is used to generate a build-specific color map for use by Auspice.
* default: `defaults/color_schemes.tsv`

### auspice_config
* type: string
* description: Path to an Auspice configuration JSON file used by `augur export`.
* default: `defaults/auspice_config.json`

### lat_longs
* type: string
* description: Path to a tab-delimited mapping of geographic scales (e.g., `location` ,`division`, etc.), geographic names (e.g., `King County`), and corresponding latitude and longitude values for the given place name. This mapping is manually curated by the Nextstrain team and updates regularly.
* default: `defaults/lat_longs.tsv`

### description
* type: string
* description: Path to a Markdown file containing a default description of each build that will be included in the build’s final Auspice JSON and appear in the build’s display in Auspice. Define a build-specific description with a path to that description file in `builds: <build_name> : description: <path_to_build_specific_description>.md`.

### clades
* type: string
* description: Path to [an Augur clade definition file](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations) where each row is a tab-delimited mapping of clade name to a gene, site (i.e., position), and alternate allele at that site for the corresponding clade.
* default: `defaults/clades.tsv`

### emerging_lineages
* type: string
* description: Path to [an Augur clade definition file](https://docs.nextstrain.org/en/latest/guides/bioinformatics/defining-clades.html#make-a-tsv-file-containing-your-clade-mutations) for emerging lineages of concern that may be a subset or variation of the lineages defined by the `clades` parameter or Pangolin lineages.
* default: `defaults/emerging_lineages.tsv`

## filter
* type: object
* description: Filters to apply to strain metadata and sequences prior to subsampling and tree inference. The workflow applies an implicit filter on the maximum collection dates later than today.

### min_length
* type: integer
* description: Minimum number of valid nucleotides (A, C, T, or G) for a genome to be included in the analysis by `augur filter --min-length`.
* default: `27000`

### exclude_where
* type: string
* description: Conditional tests of metadata columns used to exclude strains from the analysis by `augur filter --exclude-where`
* default: `"division='USA'"`

### exclude_ambiguous_dates_by
* type: string
* description: Level date ambiguity used to exclude strains from the analysis by `augur filter --exclude-ambiguous-dates-by`
* default: `any`
* examples:
	* `any`
	* `day`
	* `month`
	* `year`

### min_date
* type: float or string
* description: Minimum collection date for strains to include in the analysis used by `augur filter --min-date`. Dates can be numeric floating point values (e.g., `2019.74`) or ISO 8601-style strings (e.g., `2019-10-01`).
* default: `2019.74`

## frequencies
### min_date
* type: float or string
* description: Earliest date to estimate frequencies for. Dates can be numeric floating point values (e.g., `2019.74`) or ISO 8601-style strings (e.g., `2019-10-01`).
* default: without value supplied, defaults to 1 year before present

### max_date
* type: float or string
* description: Earliest date to estimate frequencies for. Dates can be numeric floating point values (e.g., `2021.5`) or ISO 8601-style strings (e.g., `2021-07-01`). Specifying `max_date` overrides `recent_days_to_censor`.
* default: without value supplied, defaults to today's date minus `recent_days_to_censor` parameter

### recent_days_to_censor
* type: integer
* description: How many days back from today's date should samples be hidden from frequencies calculations? This is in place to help with sampling bias where some regions have faster sequencing turnarounds than other regions.
* default: without value supplied, defaults to `0`

### pivot_interval
* type: integer
* description: Number of units between frequency estimates based on the units defined in the `pivot_interval_units` parameter. A “pivot” corresponds to a time point when frequencies are estimated.
* default: `1`

### pivot_interval_units
* type: string
* description: Unit of pivot interval spacing for frequency estimation.
* default: `weeks`
* examples:
	* `weeks`
	* `months`

### narrow_bandwidth
* type: float
* description: Variance of the KDE normal distribution in numeric floating point years (e.g., one month ~= 30 days ~= 0.08 years). This bandwidth value controls the smoothing of frequency estimates with higher values producing smoother estimates.
* default: `0.05`

### proportion_wide
* type: float
* description: Proportion of a second KDE normal distribution to add to each initial normal distribution already parameterized by the `narrow_bandwidth` parameter.
* default: `0.0`

### minimal_frequency
* Unused

### stiffness
* Unused

### inertia
* Unused

## genes
* type: array
* description: A list of genes for which `nextalign` should generate amino acid sequences during the alignment process. Gene names must match the names provided in the gene map from the `annotation` parameter.
* default: `["ORF1a", "ORF1b", "S", "ORF3a", "M", "N"]`

## include_hcov19_prefix
* type: boolean
* description: Prepend strain names with `hCoV-19/` per GISAID requirements for web display
* default: `false`

## inputs
* type: array
* description: A list of named input datasets to use for the workflow. Input order determines the precedence of genome sequences and metadata such that later datasets override earlier datasets. Each input must define a `name`, a path to `metadata`, and a path to sequences at one of many possible starting points. The workflow merged all input metadata and sequences into a single metadata and sequences file prior to subsampling.
* required
	* `name`
	* `metadata`
	* `sequences` or `aligned` or `masked` or `filtered`
* examples:
```yaml
inputs:
  - name: example-data
    metadata: data/example_metadata.tsv.xz
    sequences: data/example_sequences.fasta.xz
  - name: prealigned-data
    metadata: data/other_metadata.tsv.xz
    aligned: data/other_aligned.fasta.xz
  - name: prealigned-and-masked-data
    metadata: data/other_metadata.tsv.xz
    masked: data/other_masked.fasta.xz
  - name: prealigned-masked-and-filtered-data
    metadata: data/other_metadata.tsv.xz
    filtered: data/other_masked.fasta.xz
```

Valid attributes for list entries in `inputs` are provided below.

### name
* type: string
* description: Name of the current dataset. Names cannot contain spaces, as they correspond to files on the file system.
* examples:
	* `example-data`
	* `gisaid`
	* `washington`
	* `north-america`

### metadata
* type: string
* description: Path to a local or remote (S3) tab-delimited metadata file supported by Augur. Metadata can be uncompressed or compressed.
* examples:
	* `data/example_metadata.tsv`
	* `data/example_metadata.tsv.xz`
	* `s3://your-bucket/metadata.tsv.gz`

### sequences
* type: string
* description: Path to a local or remote (S3) FASTA file with **_un_aligned, _un_masked, and _un_filtered** genome sequences. Sequences can be uncompressed or compressed.
* examples:
	* `data/example_sequences.fasta`
	* `data/example_sequences.fasta.xz`
	* `s3://your-bucket/sequences.fasta.gz`

### aligned
* type: string
* description: Path to a local or remote (S3) FASTA file with **aligned, _un_masked, and _un_filtered** genome sequences. Sequences can be uncompressed or compressed.
* examples:
	* `data/aligned.fasta`
	* `data/aligned.fasta.xz`
	* `s3://your-bucket/aligned.fasta.gz`

### masked
* type: string
* description: Path to a local or remote (S3) FASTA file with **aligned, masked, and _un_filtered** genome sequences. Sequences can be uncompressed or compressed.
* examples:
	* `data/masked.fasta`
	* `data/masked.fasta.xz`
	* `s3://your-bucket/masked.fasta.gz`

### filtered
* type: string
* description: Path to a local or remote (S3) FASTA file with **aligned, masked, and filtered** genome sequences. Sequences can be uncompressed or compressed.
* examples:
	* `data/filtered.fasta`
	* `data/filtered.fasta.xz`
	* `s3://your-bucket/filtered.fasta.gz`

## localrules
* type: string
* description: Path to a Snakemake file to include in the workflow. This parameter is redundant with `custom_rules` and may be deprecated soon.

## logistic_growth
* type: object
* description: Parameters for estimation of logistic clade growth based on logit-transformed clade frequencies.
* required:
	* `delta_pivots`
	* `min_tips`
	* `min_frequency`
	* `max_frequency`

### delta_pivots
* type: integer
* description: Calculate logistic growth over the last N pivots which corresponds to N times the amount of time represented by the `pivot_interval_units` in the frequencies configuration.
* default: `6`

### min_tips
* type: integer
* description: The minimum number of tips a clade must have before its logistic growth is calculated.
* default: `50`

### min_frequency
* type: float
* description: The minimum current frequency for a clade to have its logistic growth calculated.
* default: `0.000001`

### max_frequency
* type: float
* description: The maximum current frequency for a clade to have its logistic growth calculated.
* default: `0.95`

## mask
* type: object
* description: Parameters for masking of invalid or problematic nucleotides in aligned sequences. In addition to the configurable parameters below, the workflow also always masks terminal gaps in the given alignment.
* required:
	* `mask_from_beginning`
	* `mask_from_end`
	* `mask_sites`

### mask_from_beginning
* type: integer
* description: Number of bases to mask from the beginning alignment.
* default: `100`

### mask_from_end
* type: integer
* description: Number of bases to mask from the end alignment.
* default: `50`

### mask_sites
* type: string
* description: Space-delimited string of 1-based genomic sites to mask
* default: `"13402 24389 24390"`

## partition_sequences
* Unused

## reference_node_name
* Unused

## refine
* type: object
* description: Parameters for inference of time trees with `augur refine`.
* required:
	* `root`
	* `clock_rate`
	* `clock_std_dev`
	* `coalescent`
	* `date_inference`
	* `divergence_unit`
	* `clock_filter_iqd`

### root
* type: string
* description: Rooting mechanism or strain name(s) whose sequences should be used to root the time tree. Only one or two (space-delimited) strain names are supported.
* default: `Wuhan/WH01/2019`
* examples:
	* `best`
	* `least-squares`
	* `min_dev`
	* `oldest`
	* `Wuhan/Hu-1/2019 Wuhan/WH01/2019`

### clock_rate
* type: float
* description: Fixed clock rate to use for time tree calculations.
* default: `0.0008`

### clock_std_dev
* type: float
* description: Standard deviation of the fixed `clock_rate` estimate.
* default: `0.0004`

### coalescent
* type: float or string
* description: Coalescent timescale in units of inverse clock rate (float), optimized as a scalar (“opt”), or skyline (“skyline”).
* default: `skyline`
* examples:
	* `opt`
	* `skyline`

### date_inference
* type: string
* description: Assign internal nodes to their jointly or marginally most likely dates.
* default: `marginal`
* examples:
	* `marginal`
	* `joint`

### divergence_unit
* type: string
* description: Units in which sequence divergence is reported.
* default: `mutations`
* examples:
	* `mutations`
	* `mutations-per-site`

### clock_filter_iqd
* type: integer
* description: Remove tips that deviate more than this number of interquartile ranges from the root-to-tip by time regression.
* default: `4`

### keep_polytomies
* type: boolean
* description: Do not attempt to resolve polytomies.
* default: `false`

### no_timetree
* type: boolean
* description: Do not produce a time tree.
* default: `false`

## run_pangolin
* type: boolean
* description: Enable annotation of Pangolin lineages for a given build’s subsampled sequences.
* default: `false`

## deploy_url
* type: string
* description: URL to an S3 bucket where Auspice JSONs should be uploaded by the `deploy` rule of the Nextstrain workflows. Only valid for Nextstrain builds.

## slack_channel
* type: string
* description: Slack channel to notify when Nextstrain builds start, fail, or get deployed. Only valid for Nextstrain builds.

## slack_token
* type: string
* description: [Slack authentication token](https://api.slack.com/authentication/token-types) required for the Slack API calls to notify the defined `slack_channel`. Only valid for Nextstrain builds.

## strip_strain_prefixes
* type: array
* description: A list of prefixes to strip from strain names in metadata and sequence records to maintain consistent strain names when analyzing data from multiple sources.
* default: `["hCoV-19/", "SARS-CoV-2/"]`

## sanitize_metadata
* type: object
* description: Parameters to configure how to sanitize metadata to a Nextstrain-compatible format.

### parse_location_field
* type: string
* description: Field in the metadata that stores GISAID-formatted location details (e.g., `North America / USA / Washington`) to be parsed into `region`, `country`, `division`, and `location` fields.
* default: `Location`

### rename_fields
* type: array
* description: List of key/value pairs mapping fields in the input metadata to rename to another value in the sanitized metadata.
* default:
```yaml
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
```

## subsampling
* type: object
* description: Schemes for subsampling data prior to phylogenetic inference to avoid sampling bias or focus an analysis on specific spatial and/or temporal scales. [See the SARS-CoV-2 tutorial for more details on defining subsampling schemes](https://docs.nextstrain.org/en/latest/tutorials/SARS-CoV-2/steps/customizing-analysis.html#subsampling).

Each named subsampling scheme supports the following attributes that the workflow passes to `augur filter`.

### group_by
* type: string
* description: Space-delimited list of metadata columns to group records by prior to subsampling to the requested or calculated number of sequences per group.
* examples:
	* `year month`
	* `region year month`

### seq_per_group
* type: integer
* description: Number of sequences to select per group of records in groups specified by `group_by`. The total number of sequences selected for each subsampling rule will be no more than the number of groups times this number of sequences per group. This parameter must be used with the `group_by` parameter.

### max_sequences
* type: integer
* description: Maximum number of sequences to select for the current subsampling rule. When used with the `group_by` parameter, Augur will calculate the number of sequences per group. When used without the `group_by` parameter, Augur will select this number of sequences at random from all available sequences. When probabilistic sampling is enabled by the `sampling_scheme` parameter, the total number of strains actually selected will be more or less than this value due to the underlying Poisson sampling process.

### sampling_scheme
* type: string
* description: A flag to pass to `augur filter` that specifies whether to enable probabilistic sampling or not. Probabilistic sampling is useful when there are more groups than requested sequences.
* default: `--probabilistic-sampling` (Augur’s default)
* examples:
	* `--probabilistic-sampling`
	* `--no-probabilistic-sampling`

### exclude
* type: string
* description: Argument to pass to `augur filter` to exclude records based on specific values in metadata columns. This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below.
* examples:
	* `"--exclude-where 'region!=Africa'"`
	* `"--exclude-where 'region!={region}'"`

### include
* type: string
* description: Argument to pass to `augur filter` to include records based on specific values in metadata columns regardless of other filters applied during subsampling (i.e., strains for which the include test evaluates to true will always be included if they exist in the metadata and sequences). This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below.
* examples:
	* `--include-where 'region=Africa'`
	* `--include-where 'region={region}'`

### query
* type: string
* description: Argument to pass to `augur filter` to select specific records by testing values in metadata columns. This argument can refer to build-specific attributes with curly bracket notation as shown in the examples below. Query values support [pandas Dataframe query syntax](https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#indexing-query) treating the metadata as a data frame.
* examples:
	* `--query "division == 'Washington'"`
	* `--query "division == '{division}'"`
	* `--query "(country == '{country}') & (division == '{division}')"`
	* `--query "division != '{division}'"`

### exclude_ambiguous_dates_by
* type: string
* description:  Level date ambiguity used to exclude strains from the analysis by `augur filter --exclude-ambiguous-dates-by`
* examples:
	* `any`
	* `day`
	* `month`
	* `year`

### min_date
* type: string
* description: Argument to `augur filter` to set the minimum collection date for strains to include in the subsampling set. Dates can be numeric floating point values (e.g., `2019.74`) or ISO 8601-style strings (e.g., `2019-10-01`).
* examples:
  * `--min-date 2019-10-01`
  * `--min-date 2019.74`

### max_date
* type: string
* description: Argument to `augur filter` to set the maximum collection date for strains to include in the subsampling set. Dates can be numeric floating point values (e.g., `2019.74`) or ISO 8601-style strings (e.g., `2019-10-01`).
* examples:
  * `--max-date 2021-04-01`
  * `--max-date 2021.25`

### priorities
* type: object
* description: Parameters to prioritize strains selected for the current subsampling rule. Currently, the workflow supports two `type`s of priority, `proximity` and `file`.
* description [proximity]: `proximity` selects samples that are genetically similar to the `focus` sample set; the `focus` sample set must be a rule in the current subsampling scheme.
* example [proximity]:
```yaml
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
```
* description [file]: `file` selects samples based on arbitrarily-defined rankings in a TSV file formatted as `strain\tnumber`. The numbers are only used to sort the samples, and are therefore arbitrary. Higher values = higher priority.

* example [file]:
```yaml
subsampling:
  my-scheme:
    my-first-rule:
      max_sequences: 10
      group_by: "country"
      priorities:
        type: "file"
        file: "path/to/priorities.tsv"
```

```
hCoV-19/USA/CZB-1234/2021	8.2
hCoV-19/USA/CZB-2345/2021	0
hCoV-19/USA/CZB-3456/2021	-3.1
```

## title
* type: string
* description: Title to provide to `augur export` and display as the title of the analysis in Auspice.

## traits
* type: object
* description: Parameters for inference of ancestral traits by `augur traits` with support for default traits and build-specific traits.
* examples:
```yaml
traits:
  default:
    sampling_bias_correction: 2.5
    columns: ["country_exposure"]
  washington:
    # Override default sampling bias correction for
    # "washington" build and continue to use default
    # trait columns.
    sampling_bias_correction: 5.0
```

Each named traits configuration (`default` or build-named) supports the following attributes.

### sampling_bias_correction
* type: float
* description: A rough estimate of how many more events would have been observed if sequences represented an even sample. [See the documentation for `augur traits` for more details](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/traits.html).
* default: `2.5`

### columns
* type: array
* description: A list of columns from the metadata for which ancestral trait values should be inferred for ancestral nodes.
* default: `["country_exposure"]`

## tree
* type: object
* description: Parameters for phylogenetic inference by `augur tree`. The tree “method” is hardcoded to `iqtree`.

### tree-builder-args
* type: string
* description: Arguments specific to the tree method (`iqtree`) to be passed through to the tree builder command run by `augur tree`.
* default: `'-ninit 10 -n 4'`

## auspice_json_prefix
* type: string
* description: Prefix to use for Auspice JSON outputs. Change this value to produce JSONs named like `auspice/<your_prefix>_global.json` for a build named `global`, for example. If you are using [Nextstrain's Community Sharing](https://docs.nextstrain.org/en/latest/guides/share/community-builds.html) to view your builds, set this value to your GitHub repository name and the `ncov` default. For example, if your repository is named `evolution`, set `auspice_json_prefix: evolution_ncov` to get JSONs you can view your `global` build at https://nextstrain.org/community/*your_github_organization*/evolution/ncov/global.
* default: `ncov`
