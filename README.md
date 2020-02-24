This is a [Nextstrain](https://nextstrain.org) build for novel coronavirus (nCoV), visible at [nextstrain.org/ncov](https://nextstrain.org/ncov).

## Data

The nCoV genomes were generously shared via [GISAID](https://gisaid.org). We gratefully acknowledge the Authors, Originating and Submitting laboratories of the genetic sequence and metadata made available through GISAID on which this research is based. For a full list of attributions please see the [metadata file](data/metadata.tsv).

nCoV genomes are not included as part of this repo as many of them are protected by the terms of GISAID sharing. These genomes will need to be supplemented by the user. Please add these as strains in `data/sequences.fasta`. Metadata for these viruses already exists in `data/metadata.tsv`. This `data/sequences.fasta` should just have strain name in the FASTA header, like so:
```
>Wuhan/WIV04/2019
attaaaggtttat...
>USA/IL1/2020
attaaaggtttat...
>Wuhan/WIV06/2019
ccttcccaggtaa...
```

## Situation Report Translations

We welcome translations of the situation reports (narratives) into languages other than English (in particular to countries affected by the outbreak), and have been very impressed with the contributions provided so far.
Please get in touch if you can help.

We suggest creating a branch for each language after the each release of the English version.
Unfortunately this means that the changes are not visible through nextstrain.org until release, but we are working on improving this.

The situation reports are generated from Markdown files (such as [this one](https://github.com/nextstrain/ncov/blob/master/narratives/ncov_sit-rep_2020-01-25.md) for 2020-01-25).

#### Current translations:

| Language | Translator(s) | Latest version released |
| -------- | ------------- | ----------------------- |
| Mandarin | [Alvin X. Han](https://twitter.com/AlvinXHan), [Fengjun Zhang](https://twitter.com/fengjun_zhang), Wei Ding | 2020-01-30 |
| Spanish  | [Ch. Julian Villabona-Arenas](https://twitter.com/Chjulian) | 2020-01-30 |
| Portugese  | [Glaucio Santos](https://twitter.com/glauciomarcos), [Anderson Brito](https://twitter.com/AndersonBrito_) | 2020-01-30 |
| French  | Etienne Simon-Lorière, Pierre Barrat-Charlaix | 2020-01-30 |
| German  | Vielen Dank, [Nicola Müller](https://twitter.com/nicfelm), [Richard Neher](https://twitter.com/richardneher) | 2020-01-30 |
| Russian  | [Ivan Aksamentov](https://twitter.com/ivan_aksamentov), Vadim Puller | 2020-01-30 |


## Contributing

We welcome contributions from the community to make this effort as useful as possible to as many people as possible.
If you spot errors or inaccuracies, please file an issue or make a pull request.
Or get in touch over email at hello@nextstrain.org or on Twitter at [@nextstrain](https://twitter.com/nextstrain).


## Building

In order to run the nextstrain build you must provision `data/sequences.fasta` and `data/metadata.tsv`.
Due to data sharing agreements, we are not able to publicaly provide these files (please see [this issue](https://github.com/nextstrain/ncov/issues/53) as we are trying to improve this).


After this, the entire build can be regenerated by running

```
snakemake -p
```
with a [local Nextstrain installation](https://nextstrain.org/docs/getting-started/local-installation) or by running
```
nextstrain build .
```
with a [containerized Nextstrain installation](https://nextstrain.org/docs/getting-started/container-installation).

The resulting output JSON at `auspice/ncov.json` can be visualized by running `auspice view --datasetDir auspice` or `nextstrain view auspice/` depending on local vs containerized installation.

_This requires [Augur](https://github.com/nextstrain/augur) version >=6.3.0, released Feb 13, 2020._


## Development

To develop locally, install required dependencies into a conda environment named `ncov` with:
```
conda env create -f environment.yml
```


## Notes

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of 5 &times; 10^-4 subs per site per year. There were SNPs present in the nCoV samples in the first and last few bases of the alignment that were masked as likely sequencing artifacts.
