# Developer guide  

## Setup  

Please see [the main Nextstrain docs](https://nextstrain.org/docs/getting-started/introduction#open-source-tools-for-the-community) for instructions for installing the Nextstrain bioinformatics pipeline (Augur) and visualization tools (Auspice).

## Data

In order to run the Nextstrain build you must provision `./data/sequences.fasta` and `./data/metadata.tsv`.
We've included a test set of sequences that are publicly available via Genbank as `./example_data/sequences.fasta`.

## Running

Please see [these docs](./docs/running.md) for instructions on how to run this build yourself.

The resulting output JSON at `auspice/ncov.json` can be visualized by running `auspice view --datasetDir auspice` or `nextstrain view auspice/` depending on local vs containerized installation.

### Finalizing automated builds

To run a regional build, be sure to update the list of regions in `nextstrain_profiles/nextstrain/builds.yaml`.
You can run all builds in parallel with the following command.

```bash
snakemake --profile nextstrain_profiles/nextstrain all_regions
```

Or you can specify final or intermediate output files like so:

```bash
# subsampled regional focal
snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_europe.json

# subsampled global
snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_global.json
```

To update ordering/lat_longs after AWS download:

```bash
snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
snakemake --profile nextstrain_profiles/nextstrain export_all_regions
```

When done adjusting lat-longs & orders, remember to run the following command to produce the final Auspice files.

```bash
snakemake --profile nextstrain_profiles/nextstrain all_regions
```
