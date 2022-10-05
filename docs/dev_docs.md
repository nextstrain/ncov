# Developer guide

## Contents

 1. [Setup](#setup)
 1. [Data](#data)
 1. [Running](#running)
 1. [Releasing new workflow versions](#releasing-new-workflow-versions)

## Running

Visit [the workflow documentation](https://docs.nextstrain.org/projects/ncov) for instructions on how to set up and run the workflow.

## Releasing new workflow versions

We use semantic versioning of the ncov workflow, denoting backward incompatible changes with major versions.
Prior to merging a pull request that introduces a new backward incompatible change (e.g., requirement of a new version of Augur), take the following steps to document these changes:

 1. Determine the new version number by incrementing [the current version](https://github.com/nextstrain/ncov/releases/) (e.g., "v2" from "v1").
 2. As part of the pull request, document the change(s) from the pull request in [`docs/src/reference/change_log.md`](https://github.com/nextstrain/ncov/blob/master/docs/src/reference/change_log.md) with the current date and new version number.
 3. Merge the pull request
 4. [Create a new GitHub release](https://github.com/nextstrain/ncov/releases/new) using the new version as the tag (e.g., "v2") and release title. Copy the changelog section for this version into the release description, along with a permalink to the changelog section (e.g. https://github.com/nextstrain/ncov/releases/v12).

We do not release new minor versions for new features, but you should document new features in the change log as part of the corresponding pull request under a heading for the date those features are merged.


## Running Core Nextstrain Builds

The "core" nextstrain builds consist of a global analysis and six regional analyses, performed independently for GISAID data and open data (currently open data is GenBank data).
Stepping back, the process can be broken into three steps:
1. Ingest and curation of raw data. This is performed by the [ncov-ingest](https://github.com/nextstrain/ncov-ingest/) repo and resulting files are uploaded to S3 buckets.
2. Phylogenetic builds, which start from the files produced by the previous step. This is performed by the profiles `nextstrain_profiles/nextstrain-open` and `nextstrain_profiles/nextstrain-gisaid`. The resulting files are uploaded to S3 buckets by the `upload` rule. 


### Manually running phylogenetic builds

To run these pipelines locally, without uploading the results:
```sh
snakemake -pf all --profile nextstrain_profiles/nextstrain-open
snakemake -pf all --profile nextstrain_profiles/nextstrain-gisaid
```
You can replace `all` with, for instance, `auspice/ncov_open_global.json` to avoid building all regions.
The resulting dataset(s) can be visualised in the browser by running `auspice view --datasetDir auspice`.

If you wish to upload the resulting information, you should run the `upload` and/or `deploy` rules.
The `upload` rule uploads the resulting files, including intermediate files, to specific S3 buckets; this rule uses the `S3_DST_BUCKET` config parameter.
The `deploy` rule uploads the dataset files such that they are accessible via nextstrain URLs (e.g. nextstrain.org/ncov/gisaid/global); this rule uses the `deploy_url` and `auspice_json_prefix` parameters.
You may wish to overwrite these parameters for your local runs to avoid overwriting data which is already present.
For instance, here are the commands used by the trial builds action (see below):
```sh
snakemake -pf upload deploy \
    --profile nextstrain_profiles/nextstrain-open \
    --config \
        S3_DST_BUCKET=nextstrain-staging/files/ncov/open/trial/TRIAL_NAME \
        deploy_url=s3://nextstrain-staging/ \
        auspice_json_prefix=ncov_open_trial_TRIAL_NAME
snakemake -pf upload deploy \
    --profile nextstrain_profiles/nextstrain-gisaid \
    --config \
        S3_DST_BUCKET=nextstrain-ncov-private/trial/TRIAL_NAME \
        deploy_url=s3://nextstrain-staging/ \
        auspice_json_prefix=ncov_gisaid_trial_TRIAL_NAME
```


## Triggering routine builds

Typically, everything’s triggered from the  `ncov-ingest` pipeline’s `trigger` command.
After updating the intermediate files, that command will run the phylogenetic `ncov` pipelines (step 3, above) force-requiring the rules `deploy` and `upload`.

## Triggering trial builds

This repository contains GitHub Actions `rebuild-gisaid` and `rebuild-open` which can be manually run [via github.com](https://github.com/nextstrain/ncov/actions).
These will run the respective phylogenetic build pipelines starting from the preprocessed (filtered) files.
This will ask for an optional “trial name” and upload intermediate files to  `nextstrain-ncov-private/trial/$TRIAL_NAME` and `nextstrain-staging/files/ncov/open/trial/$TRIAL_NAME`; if you don't supply this you will overwrite the files at `nextstrain-ncov-private` and `nextstrain-data/files/ncov/open`, as well as the trees at `nextstrain.org/ncov/gisaid/REGION` and `nextstrain.org/ncov/open/REGION`
The GitHub action will follow along with the AWS job so that you can monitor the progress; as of October 2021 each action took around 3 hours.

If you want to test a particular branch, you can select the branch to use for the trial build when running the Github action.