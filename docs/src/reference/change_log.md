# Change Log

As of April 2021, we use major version numbers (e.g. v2) to reflect backward incompatible changes to the workflow that likely require you to update your Nextstrain installation.
We also use this change log to document new features that maintain backward compatibility, indicating these features by the date they were added.

## New features since last version update

- October 21: Implement RGB-level coloring for BA.2 (21L) descendants. For background on this and lineage definitions please see [Variant report 2022-10-17](https://github.com/neherlab/SARS-CoV-2_variant-reports/blob/main/reports/variant_report_latest_draft.md). [PR 1018](https://github.com/nextstrain/ncov/pull/1018).

## v12 (12 July 2022)

- 1 July 2022: Update workflow to support Nextclade v2 (PRs [963](https://github.com/nextstrain/ncov/pull/963), [969](https://github.com/nextstrain/ncov/pull/969)). See [the Nextclade changelog](https://github.com/nextstrain/nextclade/blob/bdfd9cff73f8181bb5891a9a9c49eb0218e7e868/CHANGELOG.md#nextclade-200) for more information.

- 1 June 2022: Add "2m" timespan in Nextstrain profile builds. [PR 957](https://github.com/nextstrain/ncov/pull/957)

- 29 April 2022: Include multiple timespans in Nextstrain profile builds. [PR 910](https://github.com/nextstrain/ncov/pull/910)

- 29 April 2022: Update default mask parameters to mask 200 bases from the end of the genome rather than the existing 50. This was necessary because there is a large deletion in this region in circulating 21L viruses. This deletion is causing problems with alignment and the resulting mis-alignment appears as excess mutations in the tree. [PR 939](https://github.com/nextstrain/ncov/pull/939).

- 27 April 2022: Include new clades 22A, 22B and 22C, where 22A corresponds to Pango lineage BA.4, 22B corresponds to Pango lineage BA.5 and 22C corresponds to Pango lineage BA.2.12.1. Please see [PR 933](https://github.com/nextstrain/ncov/pull/933) for rationale behind these clade updates.

- 27 April 2022: Convert to hierarchical clade definitions. This streamlines clade definitions significantly and makes it easier to understand clade relationships. Changes can be seen in `defaults/clades.tsv` and in [PR 855](https://github.com/nextstrain/ncov/pull/855). **This feature requires Augur v14.0 or above.** To upgrade Augur follow the installation guide at [docs.nextstrain.org](https://docs.nextstrain.org/en/latest/install.html).

- 12 April 2022: Add support for numbers in build names. [PR 524](https://github.com/nextstrain/ncov/pull/524)

- 11 April 2022: Update clade definitions to be robust to presence of lineage BA.4 and BA.5 viruses. [PR #913](https://github.com/nextstrain/ncov/pull/913)

- 17 March 2022: Add `Nextclade_pango` column to metadata [PR 892](https://github.com/nextstrain/ncov/pull/892)

- 8 Mar 2022: Support disabling clock filters in the refine step by setting `clock_filter_iqd: 0` in the `refine` section. [PR #884](https://github.com/nextstrain/ncov/pull/884), [Issue #852](https://github.com/nextstrain/ncov/issues/852).

- 11 February 2022: Add colors to default Auspice config for Nextclade quality control columns and a filter for overall Nextclade QC status. [PR #861](https://github.com/nextstrain/ncov/pull/861).

## v11 (3 February 2022)

- 3 February 2022: Run Nextclade QC and filtering on the final sample set before building a tree. Nextclade also runs `nextalign` under the hood. Importantly, this enables filtering the final sample set to omit strains with many reversions and/or possible contaminants, significantly improving the quality of Omicron trees. [See the original pull request for more details](https://github.com/nextstrain/ncov/pull/842). To disable this filtering by Nextclade quality control metrics, set `skip_diagnostics: true` in [the `filter` section of your build configuration file](https://docs.nextstrain.org/projects/ncov/en/latest/reference/configuration.html#filter).

- 29 January 2022: Update "mutational fitness" coloring based on latest results from [Obermeyer et al model](https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1) via [github.com/broadinstitute/pyro-cov/](https://github.com/broadinstitute/pyro-cov/blob/master/paper/mutations.tsv).

- 26 January 2022: Template build descriptions in export for Auspice by substituting `${BUILD}` in the Markdown with the current build name.

## v10 (5 January 2022)

- 5 January 2022: Move filter and diagnostic steps after subsampling. For workflows with subsampling that does not depend on priority calculations, these changes allow the workflow to start subsampling from the metadata, skipping sequence alignment of the full input sequences and only looping through these input sequences once per build when subsampled sequences are extracted. To skip the alignment step, define your input sequences with the `aligned` directive. If you use priority-based subsampling, define your input sequences with the `sequences` directive. This reorganization of the workflow causes a breaking change in that the workflow no longer supports input-specific filtering with the `exclude_where`, `min_date`, and `exclude_ambiguous_dates_by` parameters. The workflow continues to support input-specific filtering by `min_length` and skipping of diagnostic filters with `skip_diagnostics`. [PR #814](https://github.com/nextstrain/ncov/pull/814).

- 20 December 2021: Surface the crowding penalty parameter via the config file: [PR #828](https://github.com/nextstrain/ncov/pull/827), [Issue #708](https://github.com/nextstrain/ncov/issues/708). The crowding penalty, used when calculating `priority scores` during subsampling, decreases the number of identical samples that are included in the tree during random subsampling to provide a broader picture of the viral diversity in your dataset. However, you may wish to set `crowding_penalty = 0.0` (default value = `0.1`) if you are interested in seeing as many samples as possible that are closely related to your `focal` set. You can change this parameter via `config['priorities']['crowding_penalty']`. There is no change to default behavior.

- 7 October 2021: Automatically exclude sequences whose dates are more than two weeks (this parameter can be adjusted) before the first known sequences for their corresponding clade: [PR #740](https://github.com/nextstrain/ncov/pull/740). To be extended to provide a tsv of all auto-excluded sequences and reasons for exclusion such as the reason addressed in the afformentioned PR, AKA `check_clade_dates`, and others checked in [scripts/diagnostic.py](https://github.com/nextstrain/ncov/blob/master/scripts/diagnostic.py).

- 6 October 2021: Add three configuration parameters to control the metadata sanitizer step of the workflow. These parameters allow users to specify the metadata columns to use for strain names (`metadata_id_columns`) and to resolve duplicate records with database ids (`database_id_columns`). The new `error_on_duplicate_strains` parameter allows users to ask the workflow to exit with an error when any duplicates appear in the metadata. [See the configuration reference for more details](https://docs.nextstrain.org/projects/ncov/en/latest/reference/configuration.html#sanitize-metadata). ([#728](https://github.com/nextstrain/ncov/pull/728))

- 6 October 2021: Update clades with `21I (Delta)` and `21J (Delta)` viruses. These are subclades within `21A (Delta)`. Based on mutations they should have largely Delta-like phenotypes, although additional ORF1a mutations in `21J (Delta)` appear to confer higher fitness.

## v9 (6 October 2021)

- 6 October 2021: Remove travel exposure adjustment. This is a potentially breaking change if you had explicitly opted into this functionality (by default, it was disabled); in this case your snakefile will print an appropriate error and exit. Note that you can still define exposure metadata traits as colorings in your `auspice-config.json`, however be aware that these may not be curated in the future. See [PR 723](https://github.com/nextstrain/ncov/pull/723) for more.

- 6 October 2021: Change how the internal (core) nextstrain profiles are run. We now split each (GISAID, open) into a preprocessing profile and a phylogenetic build pipeline. Please see [PR 730](https://github.com/nextstrain/ncov/pull/730) for more. There should be no changes to other profiles.

- 19 September 2021: Include "mutational fitness" coloring based on [Obermeyer et al model](https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1). This annotates each node in the tree with a `mutational_fitness` trait by summing mutational effects from Obermeyer et al Supplementary Data S2. This should provide a mechanism to flag emergence of novel variants that may have higher fitness than circulating viruses.

## v8 (19 Aug 2021)

### Major changes

- Annotate CDC-style epiweeks (e.g., "202019") as a color-by and filter option in Auspice JSONs ([#703](https://github.com/nextstrain/ncov/pull/703)). This functionality requires [the Python epiweeks package](https://pypi.org/project/epiweeks/). You will need to update your software environment to include this package, depending on how you run your builds.
  - If you use the Nextstrain CLI with Docker, update the Docker image with `nextstrain update` and then run your builds as usual with `nextstrain build`.
  - If you use the Nextstrain CLI without Docker, run your builds with `nextstrain build . --use-conda <...other options...>`.
  - If you use Snakemake, run your builds with `snakemake --use-conda <...other options...>`.
  - If you manage your own Conda environment, install epiweeks manually in the environment with `conda install -c bioconda epiweeks`.

### Features

- Update Conda environment to use [Augur 13.0.0](https://github.com/nextstrain/augur/blob/master/CHANGES.md#1300-17-august-2021) for an improved filtering experience ([#703](https://github.com/nextstrain/ncov/pull/703)).

## New features since last version update

 - 11 August 2021: Add support for "Sequences" and "Patient status metadata" downloads from GISAID's search interface including [documentation in the tutorial of how to use these data](../guides/data-prep.html#curate-data-from-gisaid-search-and-downloads). ([#701](https://github.com/nextstrain/ncov/pull/701))
 - 6 August 2021: We've replaced the mechanisms that support remote file inputs (e.g. `s3://` URLs) to improve internal workflow structure, extend support to `gs://`, `http://`, and `https://` URLs, and expand support for compressed inputs.
   Our [remote file inputs documentation](remote_inputs) is updated to reflect the changes.

   This change should be backwards compatible and largely transparent to end users.
   The most visible change for anyone using remote file inputs is the local download location of the remote files: instead of being within the `results/` directory, dynamic directories based on the remote URL are now used.
   For example, the remote `metadata` input `https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz` would be downloaded to the local path `data.nextstrain.org/files/ncov/open/metadata.tsv.gz`; the remote `sequences` input `s3://nextstrain-data/files/ncov/open/sequences.fasta.xz` would be downloaded to `nextstrain-data/files/ncov/open/sequences.fasta.xz`.

 - 7 July 2021: Extensive changes to internally-used Nextstrain profiles.
 There should be no breaking changes to existing workflows outside of `./nextstrain_profiles/`.
 Please see [PR #628](https://github.com/nextstrain/ncov/pull/628) for full details. Briefly:
     - The (GISAID) profile has been renamed to `./nextstrain_profiles/nextstrain-gisaid`
     - A new "open" (GenBank) profile has been added `./nextstrain_profiles/nextstrain-open`
     - Intermediate open (GenBank) files, including sequences, & alignments are now publicly available for workflows to use as starting points. See the [remote inputs documentation](remote_inputs) for details.
 - 3 July 2021: Allow optional prefixing of `hCoV-19/` to strain names when exporting Auspice JSON for visualization. This is specified via the config option `include_hcov19_prefix`. This is included in Nextstrain-maintained builds at the request of GISAID.
 - 27 June 2021: Update clade definitions with 21G (Lambda, C.37) emerging from Peru and 21H (B.1.621) emerging from Colombia.
 - 22 June 2021: Add the ability to specify subsampling via a `priorities.tsv` file. To use, set the `priorities > type: file` and add `priorities > file: path/to/priorities.tsv` to your build's `subsampling` schema. `priorities.tsv` contains `strain name\tarbitrary numerical value`. Higher values = higher priority. ([#664](https://github.com/nextstrain/ncov/pull/664))
 - 18 June 2021: Change default behavior of frequency estimation to estimate frequencies starting 1 year prior to the current date. To override this default behavior, define a `min_date` in the `frequencies` section of the builds configuration. ([#659](https://github.com/nextstrain/ncov/pull/659))

## v7 (27 May 2021)

For more details about this release, see [the configuration reference for the new "sanitize metadata" parameters](configuration.html#sanitize_metadata) and [the corresponding pull request](https://github.com/nextstrain/ncov/pull/640).

### Major changes

- Deduplicate metadata and sequences from each `inputs` dataset at the beginning of the workflow.

### Features

- Support full GISAID metadata and sequences from the "Download packages" interface by converting this default format into Nextstrain-compatible metadata and sequences.
- Support reading metadata and sequences directly from GISAID's tar archives. For example, you can now define `inputs` as `metadata: data/ncov_north-america.tar.gz` and `sequences: data/ncov_north-america.tar.gz` to decompress and read the corresponding data from the archive.

## New features since last version update

 - 25 May 2021: Support custom Auspice JSON prefixes with a new configuration parameter, `auspice_json_prefix`. [See the configuration reference for more details](configuration.html#auspice_json_prefix). ([#643](https://github.com/nextstrain/ncov/pull/643))

## v6 (20 May 2021)

### Major changes

- Fix bug in precedence of input data such that duplicate sequence and metadata records are resolved by always preferring the record from the last `inputs` dataset. Thank you to @ttung for catching/patching this! If you have depended on the previous behavior where the sequence from first `input` dataset was preferred, you will need to change the order of your `inputs` such that the preferred input appears last in the list. ([#639](https://github.com/nextstrain/ncov/pull/639)).

## New features since last version update

 - 19 May 2021: Compress metadata, sequence indices, and early intermediate sequences (aligned, masked, filtered, combined for subsampling, and subsampled files) to save disk space. ([#636](https://github.com/nextstrain/ncov/pull/636))
 - 12 May 2021: Include S1 mutations and nextalign-based ancestral amino acid mutations in Auspice JSONs by default instead of requiring the now-unnecessary `use_nextalign` configuration parameter. ([#630](https://github.com/nextstrain/ncov/pull/630))
 - 12 May 2021: [Document all available workflow configuration parameters](https://nextstrain.github.io/ncov/configuration). ([#633](https://github.com/nextstrain/ncov/pull/633))

## v5 (7 May 2021)

[See the corresponding pull request](https://github.com/nextstrain/ncov/pull/615) for more details about this release.

### Major changes

- Drop support for old sequence/metadata inputs. This change removes support for the `config["sequences"]` and `config["metadata"]` starting points for the workflow in favor of the more flexible [`config["inputs"]` format](configuration.html#inputs).
- Use `nextalign` for alignment instead of `mafft`. This change completely removes support for `mafft` in favor of `nextalign`. Future versions may reinstate `mafft` support as part of `augur align` updates.

### Minor changes

- Drop unused haplotype status rule and script
- Remove unused nucleotide mutation frequencies rule
- Use `augur distance` for mutation counts instead of a custom script in the ncov repository. [Recent improvements to `augur distance` in v12.0.0](https://github.com/nextstrain/augur/blob/master/CHANGES.md#1200-13-april-2021) enable this change by properly accounting for insertion/deletion events.

## v4 (5 May 2021)

[See the corresponding pull request](https://github.com/nextstrain/ncov/pull/605) for more details about changes in this release.

### Major changes

- Change the default build name from "global" to "default-build" and use a default subsampling scheme that selects all input sequences
- Warn about duplicate sequences found when merging sequences from multiple inputs instead of throwing an error (set `combine_sequences_for_subsampling: warn_about_duplicates: false` in your configuration file to revert this behavior)

### Features

- Define a new subsampling scheme named `all` that selects all input sequences
- Add a new top-level configuration parameter `default_build_name` to allow overriding new default name of "default-build"
- Support compressed sequence inputs for alignment with mafft and nextalign (requires mafft upgrade)
- Sanitize strain names in sequences and metadata from different sources (e.g., `hCoV-19/` from GISAID or `SARS-CoV-2/` from GenBank, etc.)

## New features since last version update

- 20 April 2021: Surface emerging lineage as a colorby. This replaces the rather stale color by "Emerging Clade" with a new color by "Emerging Lineage". This focuses on PANGO lineages that are of interest triangulated by [CoVariants](https://covariants.org/), [PANGO](https://cov-lineages.org/) international lineage reports, [CDC](https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/variant-surveillance/variant-info.html) VUIs and VOCs and [PHE](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/975742/Variants_of_Concern_VOC_Technical_Briefing_8_England.pdf) VUIs and VOCs. The intention is for the listing at `emerging_lineages.tsv` to be updated frequently with new lineages added and no longer interesting lineages dropped. [#609](https://github.com/nextstrain/ncov/pull/609)

- 12 April 2021: Calculate current clade frequency and logistic growth rate across nodes in the phylogeny. This produces a new `logistic_growth.json` file and uses this file to add a coloring the final Auspice JSON. Implementation choices are discussed in PR [#595](https://github.com/nextstrain/ncov/pull/595).

- 12 April 2021: Annotate Pangolin lineages per build in a `pangolineages.json` file and final Auspice JSON by adding `run_pangolin: true` to the top-level workflow config (`builds.yaml`). **Note: this annotation only works when running the workflow with Snakemake's `--use-conda` flag or if your environment has [Pangolin installed](https://github.com/cov-lineages/pangolin).** [#593](https://github.com/nextstrain/ncov/pull/593)


## v3 (12 April 2021)

- Use Augur 11.2.0's metadata-only output arguments to aggregate subsampled sequences and metadata [#592][]
- Use Augur 11.3.0's `io.py` module to combine and deduplicate uncompressed or compressed sequences when handling multiple input datasets [#592][]

[#592]: https://github.com/nextstrain/ncov/pull/592

## v2 (9 April 2021)

This release reflects the state of the workflow when we instituted our workflow versioning system.

## July 2020 Update

In order to make our repository more accessible for first-time users, and those who want to create their own customized build, we made a number of changes in July 2020.
Here, we outline those changes, and also list 'breaking changes' that may effect advanced users who were already running their own builds.

You can check the lists below to see if there's any action you could take to make updating to the latest version of this repository is smooth!

### Summary of changes

The goal of this release is to make it easier for new users to get their own `ncov` build up and running. This consists of a new, extensive tutorial hosted on github pages; a simplified repository structure; and more didactic file names.

Summary of breaking changes:
- All nextstrain-specific files now live under `nextstrain_profiles`
- `config` is now named `defaults`
- `envs`, `rules` and `schemas` now live under `workflow`

#### Additions and updates
- README is updated
- Adds in-depth tutorial to `docs/` and accompanying github pages setup files
- `data/` now contains example data that should remain packaged with the repository
- `narratives/` now contains only a `template_narrative.md` for new users (not to worry -- all finished Situation Reports are still hosted over in the `nextstrain/narratives` repo)
- `my_profiles` now contains extensive `example` (formerly `king-county`) and `example_advanced_customization` (formerly `swiss`) profiles

#### File name changes
- The main `config.yaml` file is now more precisely labeled as `parameters.yaml`
- `reference.gb` > `reference_seq.gb`
- `ordering.tsv` > `color_ordering.tsv`
- `config` > `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`
- nextstrain-specific profiles labeled with `nextstrain_*`


#### Structural changes
To improve repo clarity and approachability for new users, the following files have been reorganized.
- Workflow-related files now live in a `workflow` directory, rather than the top level
  - `envs` > `workflow/envs`
  - `schemas` > `workflow/schemas`
  - `rules` > `workflow/snakemake_rules`
- Nextstrain-specific profiles now live under `nextstrain_profiles`

(Reiterated from last section because they're important :)
- `config` is now `defaults`, with nextstrain-specific files (e.g., `description.md`) moved to `nextstrain_profiles/nextstrain/`
- `config.yaml` is now `parameters.yaml`

### Breaking Changes

If you've got custom builds running and regularly pull from the repository, you should check here to see if any of the changes you've made might make the merge harder.
Usually, just copying and renaming folders and/or files is enough to avoid too many merge conflicts!

First, you should check the file locations that have changed, above.
If you reference any of these files specifically in your pipelines or profile, be sure to update those links!

Be sure to also check any links to files that were previously in `profiles`.
For example, if you were referencing the cluster information in `profiles/nextstrain-scicore/cluster.json` for example) - we've now got two profiles folders - you'd find `nextstrain-scicore` in the `nextstrain_profiles` folder now.

If you modified the `swiss`, `default`, or `king-county` profiles, this may cause a merge conflict, as these folders have been renamed.
To keep your builds as they are, we recommend copying these folders, and then after the merge, restoring the name you prefer.

We now have a `data` folder which comes with the repository and contains example build data.
If you had previously created a `data` folder, this may cause a merge conflict - usually this can be easily resolved.


### Developer Changes

For those who help run and maintain the Nextstrain.org builds, there's a few changes that might be different from what you're used to.

Just focusing on the files mostly accessed by the 'build shepherds,' here's what's good to know:
- The `config` folder is now `defaults`
- `orderings.tsv` is now `color_orderings.tsv`
- `rules/nextstrain_exports.smk` (the rules that do the 'finalizing' of the builds) isnow `workflow/snakemake_rules/export_for_nextstrain.smk`

The 'finishing builds' commands are also slightly different because the profiles have moved - they'll look familiar though!
You can already find these changed commands in the `export_for_nextstrain.smk` file, but here they are as well:

```
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
#   snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
# to produce the final Auspice files!
```
