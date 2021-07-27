# Overview of remote nCoV files (intermediate build assets)

This page provides an overview of intermediate files which Nextstrain produces.
Where appropriate, these files can be starting points for the [ncov pipeline](https://github.com/nextstrain/ncov/) (discussed below).

We have two GitHub repositories which routinely upload files to [S3 buckets](https://aws.amazon.com/s3/): [ncov-ingest](https://github.com/nextstrain/ncov-ingest/) and [ncov](https://github.com/nextstrain/ncov/).
Each of those runs separate pipelines for GISAID and GenBank (aka "open") data sources; these pipelines start with data curation and QC steps and end with the phylogenetic analyses you can see on [nextstrain.org](https://nextstrain.org/sars-cov-2/)

The GISAID data is stored at `s3://nextstrain-ncov-private` and is not publicly available, in line with the GISAID Terms of Use (this is used internally by Nextstrain).

The open (GenBank) data is publicly available at three endpoints:

  - `https://data.nextstrain.org/files/ncov/open/`
  - `s3://nextstrain-data/files/ncov/open/`
  - `gs://nextstrain-data/files/ncov/open/` (mirrored daily from S3 by the Broad Institute)

**Our intention is to make GenBank intermediate files open and available for everyone to use, and to keep these files up-to-date.**
The paths for specific files are the same under each endpoint, e.g. `https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz`, `s3://nextstrain-data/files/ncov/open/metadata.tsv.gz`, and `gs://nextstrain-data/files/ncov/open/metadata.tsv.gz` all exist.
See below for a list of files that exist.

## All available genomes and metadata
Entire metadata & sequences data is uploaded from the `ncov-ingest` workflows for each of the `gisaid` and `open` sources:

* `metadata.tsv.gz`
* `sequences.fasta.xz`
* `nextclade.tsv.gz`
* `additional_info.tsv.gz` (GISAID only)
* `flagged_metadata.txt.gz` (GISAID only)

The `ncov` repository contains two Nextstrain-only workflows which use the above data to generate our analyses visible on nextstrain.org.
A side-effect of this is the creation and upload of processed versions of the entire dataset:

* `aligned.fasta.xz` alignment via [nextalign](https://github.com/nextstrain/nextclade/tree/master/packages/nextalign_cli). The default reference genome is [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) (Wuhan-Hu-1).
* `mutation-summary.tsv.xz` A summary of the data in `aligned.fasta.xz`.
* `masked.fasta.xz` Masked alignment
* `filtered.fasta.xz` The masked alignment excluding data with incomplete / invalid dates, unexpected genome lengths, missing metadata etc. We also maintain a [list of sequences to exclude](https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt) which are removed at this step. These sequences represent duplicates, outliers in terms of divergence or sequences with faulty metadata.

## Subsampled datasets

Our GISAID and GenBank (open) profiles each define 7 builds (a Global build and one build per region: Africa, Asia, Europe, Oceania, North and South America).
Each of these is a different subsample of the entire dataset, and each will result in the following intermediates uploaded:

* `{build_name}/sequences.fasta.xz`
* `{build_name}/metadata.tsv.xz`
* `{build_name}/aligned.fasta.xz`
* `{build_name}/{build_name}.json` (the main Auspice dataset file)
* `{build_name}/{build_name}_tip-frequencies.json`
* `{build_name}/{build_name}_root-sequence.json`

---

# Summary of available GenBank (open) files

Each regional build (`global`, `africa`, `asia`, `europe`, `north-america`, `oceania` and `south-america`) contains a subsampled set of approximately 4000 sequences.
They are a good starting point if you are seeking a representative sample of data.
Where available, this table also provides the URL for the resulting Auspice visualisation of the data.

> Please note that these files are uploaded in two batches (see above for details).
This means that the full GenBank metadata and sequences are typically updated a couple of hours before the more processed files.

| description          | type      | address                                                         |
| ---                  | ---       | ---                                                             |
| Full GenBank data    | metadata  | s3://nextstrain-data/files/ncov/open/metadata.tsv.gz            |
|                      | sequences | s3://nextstrain-data/files/ncov/open/sequences.fasta.xz         |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/aligned.fasta.xz           |
|                      | masked    | s3://nextstrain-data/files/ncov/open/masked.fasta.xz            |
|                      | filtered  | s3://nextstrain-data/files/ncov/open/filtered.fasta.xz          |
| Global sample        | metadata  | s3://nextstrain-data/files/ncov/open/global/metadata.tsv.xz     |
|                      | sequences | s3://nextstrain-data/files/ncov/open/global/sequences.fasta.xz  |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/global/aligned.fasta.xz    |
|                      | auspice   | nextstrain.org/ncov/open/global                                 |
| Africa sample        | metadata  | s3://nextstrain-data/files/ncov/open/africa/metadata.tsv.xz     |
|                      | sequences | s3://nextstrain-data/files/ncov/open/africa/sequences.fasta.xz  |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/africa/aligned.fasta.xz    |
|                      | auspice   | nextstrain.org/ncov/open/africa                                 |
| Asia sample          | metadata  | s3://nextstrain-data/files/ncov/open/asia/metadata.tsv.xz       |
|                      | sequences | s3://nextstrain-data/files/ncov/open/asia/sequences.fasta.xz    |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/asia/aligned.fasta.xz      |
|                      | auspice   | nextstrain.org/ncov/open/asia                                   |
| Europe sample        | metadata  | s3://nextstrain-data/files/ncov/open/europe/metadata.tsv.xz     |
|                      | sequences | s3://nextstrain-data/files/ncov/open/europe/sequences.fasta.xz  |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/europe/aligned.fasta.xz    |
|                      | auspice   | nextstrain.org/ncov/open/europe                                 |
| North America sample | metadata  | s3://nextstrain-data/files/ncov/open/north-america/metadata.tsv.xz    |
|                      | sequences | s3://nextstrain-data/files/ncov/open/north-america/sequences.fasta.xz |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/north-america/aligned.fasta.xz   |
|                      | auspice   | nextstrain.org/ncov/open/north-america                          |
| Oceania sample       | metadata  | s3://nextstrain-data/files/ncov/open/oceania/metadata.tsv.xz    |
|                      | sequences | s3://nextstrain-data/files/ncov/open/oceania/sequences.fasta.xz |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/oceania/aligned.fasta.xz   |
|                      | auspice   | nextstrain.org/ncov/open/oceania                                |
| South America sample | metadata  | s3://nextstrain-data/files/ncov/open/south-america/metadata.tsv.xz    |
|                      | sequences | s3://nextstrain-data/files/ncov/open/south-america/sequences.fasta.xz |
|                      | aligned   | s3://nextstrain-data/files/ncov/open/south-america/aligned.fasta.xz   |
|                      | auspice   | nextstrain.org/ncov/open/south-america                          |


---

# Starting your build from these intermediates

Each workflow defines one or more inputs in the `builds.yaml` file.

In the simplest form, an input specifies a local path to some metadata and sequences, like so:

```yaml
inputs:
  - name: example-data
    metadata: data/example_metadata.tsv
    sequences: data/example_sequences.fasta
```

Using the above table, we can easily modify this to create a build which uses the global subsample of GenBank data:

```yaml
inputs:
  - name: global-representative-genbank-sample
    metadata: s3://nextstrain-data/files/ncov/open/global/metadata.tsv.gz
    sequences: s3://nextstrain-data/files/ncov/open/global/sequences.fasta.gz
```

To avoid unnecessarily aligning these sequences, we can instead start from the aligned sequences, like so:

```yaml
inputs:
  - name: global-representative-genbank-sample
    metadata: s3://nextstrain-data/files/ncov/open/global/metadata.tsv.gz
    aligned: s3://nextstrain-data/files/ncov/open/global/aligned.fasta.gz
```

The following starting points are available:

* replace `sequences` with `aligned` (skips alignment 😉)
* replace `sequences` with `masked` (skips alignment and masking steps)
* replace `sequences` with `filtered` (skips alignment, masking and basic filtering steps)


## Compressed vs uncompressed starting points

In general, we try to transparently support compressed metadata and sequences for any input stage.
Files may be compressed using `xz` (`.xz`) or `gzip` (`.gz`) compression.
The following table summarises the current situation:

|            | local (uncompressed)  | local (compressed)| remote (uncompressed) | remote (compressed)|
| ---        | ---                   | ---               | ---                   | ---               |
| sequences  |  ✅                   |  ✅                |  ❌  <sup>1</sup>     |  ✅  <sup>2</sup> |
| metadata   |  ✅                   |  ✅                |  ✅                   |  ✅  <sup>3</sup> |
| aligned    |  ✅                   |  ✅                |  ✅                   |  ✅  <sup>3</sup> |
| masked     |  ✅                   |  ✅                |  ✅                   |  ✅  <sup>3</sup> |
| filtered   |  ✅                   |  ✅                |  ✅                   |  ✅  <sup>3</sup> |

1. Download succeeds but `.xz` suffix is appended which causes errors downstream in the pipeline.
2. File is downloaded with a `.xz` suffix, but gzip compression will still work!
3. Decompressed during download.
