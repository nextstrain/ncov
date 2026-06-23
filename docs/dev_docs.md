# Developer guide

## Contents

 1. [Setup](#setup)
 2. [Data](#data)
 3. [Running](#running)
 4. [Releasing new workflow versions](#releasing-new-workflow-versions)
 5. [Proposing new clade designations](#proposing-new-clade-designations)
 6. [Running Core Nextstrain Builds](#running-core-nextstrain-builds)

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


## Proposing new clade designations

Designating a new Nextstrain clade is a judgement call on top of a quantitative threshold:
a candidate Pango lineage should carry **≥1 spike mutation** relative to its parent clade and
have risen to **>30% regional** or **>20% global** frequency (see the
[clade-naming blog post](https://nextstrain.org/blog/2022-04-29-SARS-CoV-2-clade-naming-2022)).
We keep a strict 1:1 mapping between a Pango lineage and a Nextstrain clade name, and prefer to
demarcate on the branch with the larger fitness differential.

`scripts/propose_clades.py` automates the *evidence gathering* for this. It runs entirely on
**open** (GenBank/INSDC) data — which is what we have since GISAID API access was lost (see the
[2026-02-24 blog post](https://nextstrain.org/blog/2026-02-24-gisaid-joint-response)) and doesn't
directly edit `defaults/`. It assembles a ranked candidate dossier; a human (or Claude Code, 
see below) makes the call and writes the edits.

### Running the script

```sh
python3 scripts/propose_clades.py
```

The script uses only the Python standard library. It expects the local open global build
(`auspice/ncov_open_global_6m.json` and `…_tip-frequencies.json`) to be present (this is 
purposeful as next steps rely on re-running build with updated `clades.tsv`), and fetches
the rest automatically (cached under `results/clade_cache/`, re-pulled with `--refresh`):

- the [Nextclade SARS-CoV-2 reference tree](https://nextstrain.org/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs) — clean lineage-defining mutations + parent clade;
- the North America and Europe open 6m builds — regional tip-frequencies;
- the [LAPIS open API](https://lapis.cov-spectrum.org/open/v2/) — exact per-lineage regional frequency;
- the [forecasts-ncov open MLR results](https://nextstrain.github.io/forecasts-ncov/) — growth advantage (fitness) and modeled frequency.

Useful flags: `--refresh` (ignore cache), `--skip-lapis` / `--skip-mlr` / `--skip-regional` for a
fast offline pass, and `--designation-year` to override the year used for suggested clade names
(defaults to the current year — the name's year is the *designation* year, not the lineage's
emergence year). Run `python3 scripts/propose_clades.py --help` for the full list.

Outputs are written to (git-ignored) `results/`:

- `results/clade_candidates.md` — the human-readable dossier;
- `results/clade_candidates.json` — the same data, structured, with paste-ready file edits.

Because open data is sparse and ~2 months lagged, the dossier shows each frequency alongside its
source (build tip-frequency, LAPIS, MLR) so the numbers can be triangulated and cross-checked
against the live resources rather than trusted blindly.

### Intended path for Claude Code

This tool is designed to be driven by [Claude Code](https://claude.com/claude-code) via the
`designate-clades` skill (`.claude/skills/designate-clades/SKILL.md`). The skill walks through the
full loop:

1. Run `scripts/propose_clades.py` and read `results/clade_candidates.md`.
2. For each flagged candidate, cross-check against the live resources (the open builds at
   `nextstrain.org/ncov/open/{global,north-america,europe}/6m`, the Cov-Spectrum lineage page,
   forecasts-ncov, and the Nextclade reference tree).
3. Decide the demarcation: among nested candidate lineages, elevate the one where the MLR growth
   advantage jumps most and frequency is rising; or decline if it is just slicing the current
   dominant clade finer.
4. Write the edits (see below) and validate them.
5. Draft the PR.

Open the skill in a Claude Code session in this repo and ask it to find or designate new clades.

### How a designation PR is made

A designation edits five files in `defaults/` (the dossier prints a paste-ready block per
candidate; trim the `clades.tsv` nuc rows to a characteristic handful — the spike-linked rows are
flagged for you):

1. `clades.tsv` — `<NAME>\tclade\t<PARENT>` plus a few `<NAME>\tnuc\t<site>\t<alt>` rows
2. `clade_display_names.yml` — `<NAME>: <NAME> (<Pango>)`
3. `clade_hierarchy.tsv` — `<NAME>\t<PARENT>\t<WHO>`
4. `clade_emergence_dates.tsv` — `<NAME>\t<YYYY-MM-01>` (the lineage's first-seen date)
5. `color_ordering.tsv` — `clade_membership\t<NAME> (<Pango>)`

Validate that the new definition parses and inherits correctly through Augur (needs a Python
environment with `augur` installed, e.g. the workflow's conda environment):

```sh
scripts/expand-clade-definitions defaults/clades.tsv | grep "^<NAME>"
```

Then open a pull request in the style of
[#1149](https://github.com/nextstrain/ncov/pull/1149),
[#1152](https://github.com/nextstrain/ncov/pull/1152),
[#1158](https://github.com/nextstrain/ncov/pull/1158), and
[#1192](https://github.com/nextstrain/ncov/pull/1192). The justification should state the new clade
name and its Pango lineage, the parent clade, the defining spike mutations, the frequency evidence
(with regions and numbers), the MLR growth advantage, and a one-line note on the demarcation
choice. Keep the dossier's exact numbers in the PR body so reviewers can trace them.

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


### Triggering routine builds

Typically, everything’s triggered from the  `ncov-ingest` pipeline’s `trigger` command.
After updating the intermediate files, that command will run the phylogenetic `ncov` pipelines (step 3, above) force-requiring the rules `deploy` and `upload`.

### Triggering trial builds

This repository contains GitHub Actions `rebuild-gisaid` and `rebuild-open` which can be manually run [via github.com](https://github.com/nextstrain/ncov/actions).
These will run the respective phylogenetic build pipelines starting from the preprocessed (filtered) files.
This will ask for an optional “trial name” and upload intermediate files to  `nextstrain-ncov-private/trial/$TRIAL_NAME` and `nextstrain-staging/files/ncov/open/trial/$TRIAL_NAME`; if you don't supply this you will overwrite the files at `nextstrain-ncov-private` and `nextstrain-data/files/ncov/open`, as well as the trees at `nextstrain.org/ncov/gisaid/REGION` and `nextstrain.org/ncov/open/REGION`
The GitHub action will follow along with the AWS job so that you can monitor the progress; as of October 2021 each action took around 3 hours.

If you want to test a particular branch, you can select the branch to use for the trial build when running the Github action.
