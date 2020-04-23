# Developer guide  

## Setup  

Please see [the main Nextstrain docs](https://nextstrain.org/docs/getting-started/introduction#open-source-tools-for-the-community) for instructions for installing the Nextstrain bioinformatics pipeline (Augur) and visualization tools (Auspice).

## Data

In order to run the Nextstrain build you must provision `./data/sequences.fasta` and `./data/metadata.tsv`.
We've included a test set of sequences that are publicly available via Genbank as `./example_data/sequences.fasta`.

## Running

Please see [these docs](./docs/running.md) for instructions on how to run this build yourself.

The resulting output JSON at `auspice/ncov.json` can be visualized by running `auspice view --datasetDir auspice` or `nextstrain view auspice/` depending on local vs containerized installation.

