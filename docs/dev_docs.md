# Developer guide

## Contents

 1. [Setup](#setup)
 1. [Data](#data)
 1. [Running](#running)
 1. [Releasing new workflow versions](#releasing-new-workflow-versions)

## Setup

Please see [the main Nextstrain docs](https://nextstrain.org/docs/getting-started/introduction#open-source-tools-for-the-community) for instructions for installing the Nextstrain bioinformatics pipeline (Augur) and visualization tools (Auspice).

## Data

In order to run the Nextstrain build you must provision `./data/sequences.fasta` and `./data/metadata.tsv`.
We've included a test set of sequences that are publicly available via Genbank as `./example_data/sequences.fasta`.

## Running

Please see [these docs](./docs/running.md) for instructions on how to run this build yourself.

The resulting output JSON at `auspice/ncov.json` can be visualized by running `auspice view --datasetDir auspice` or `nextstrain view auspice/` depending on local vs containerized installation.

### Finalizing automated builds

To run a regional build, be sure to update the list of regions in `nextstrain_profiles/nextstrain-gisaid/builds.yaml`.
You can run all builds in parallel with the following command.

```bash
snakemake --profile nextstrain_profiles/nextstrain-gisaid all_regions
```

Or you can specify final or intermediate output files like so:

```bash
# subsampled regional focal
snakemake --profile nextstrain_profiles/nextstrain-gisaid auspice/ncov_europe.json

# subsampled global
snakemake --profile nextstrain_profiles/nextstrain-gisaid auspice/ncov_global.json
```

To update ordering/lat_longs after AWS download:

```bash
snakemake --touch --forceall --profile nextstrain_profiles/nextstrain-gisaid
snakemake --profile nextstrain_profiles/nextstrain-gisaid clean_export_regions
snakemake --profile nextstrain_profiles/nextstrain-gisaid export_all_regions
```

When done adjusting lat-longs & orders, remember to run the following command to produce the final Auspice files.

```bash
snakemake --profile nextstrain_profiles/nextstrain-gisaid all_regions
```

## Releasing new workflow versions

We use semantic versioning of the ncov workflow, denoting backward incompatible changes with major versions.
Prior to merging a pull request that introduces a new backward incompatible change (e.g., requirement of a new version of Augur), take the following steps to document these changes:

 1. Determine the new version number by incrementing [the current version](https://github.com/nextstrain/ncov/releases/) (e.g., "v2" from "v1").
 2. As part of the pull request, document the change(s) from the pull request in [`docs/change_log.md`](https://github.com/nextstrain/ncov/blob/master/docs/change_log.md) with the current date and new version number.
 3. Merge the pull request
 4. [Create a new GitHub release](https://github.com/nextstrain/ncov/releases/new) using the new version as the tag (e.g., "v2") and release title. Leave the release description empty.

We do not release new minor versions for new features, but you should document new features in the change log as part of the corresponding pull request under a heading for the date those features are merged.


## Triggering routine builds

Typically, everything’s triggered from the  `ncov-ingest` pipeline’s `trigger` command.
After updating the intermediate files, that command will run this `ncov` pipeline force-requiring the rules `deploy` and `upload`.

## Triggering trial builds

This repository contains a GitHub Action `trial-build` which is manually run [via github.com](https://github.com/nextstrain/ncov/actions/workflows/trial-build.yml).
This will ask for a “trial name” and upload intermediate files to  `nextstrain-ncov-private/trial/$TRIAL_NAME` and `nextstrain-staging/files/ncov/open/trial/$TRIAL_NAME`.
Auspice datasets for visualisation will be available at `https://nextstrain.org/staging/ncov/gisaid/trial/$TRIAL_NAME/$BUILD_NAME` and `https://nextstrain.org/staging/ncov/open/trial/$TRIAL_NAME/$BUILD_NAME`.
