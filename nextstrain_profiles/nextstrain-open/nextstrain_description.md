Compiled Nextstrain SARS-CoV-2 resources are available at [nextstrain.org/sars-cov-2](https://nextstrain.org/sars-cov-2/). Follow [@nextstrain](https://twitter.com/nextstrain) for updates.

This phylogeny shows evolutionary relationships of SARS-CoV-2 viruses from the ongoing COVID-19 pandemic. Although the genetic relationships among sampled viruses are generally quite clear, there is considerable uncertainty surrounding estimates of specific transmission dates and in reconstruction of geographic spread. Please be aware that specific inferred geographic transmission patterns and temporal estimates are only a hypothesis.

There are millions of complete SARS-CoV-2 genomes available on open databases and this number increases every day. This visualization can only handle ~4000 genomes in a single view for performance and legibility reasons. Because of this we subsample available genome data for our analysis views. We provision multiple views to focus subsampling on different geographic regions and different time periods. These views are available through the "Dataset" dropdown on the left or by clicking on the following links:

| &nbsp;            | past 2 months                                                            | past 6 months                                                            | all time                                                                             |
| ----------------- | ------------------------------------------------------------------------ | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ |
| **global**        | [global/2m](/ncov/open/global/2m)                                        | [global/6m](/ncov/open/global/6m)                                        | [global/all-time](/ncov/open/global/all-time)                                        |
| **Africa**        | [africa/2m](/ncov/open/africa/2m?f_region=Africa)                        | [africa/6m](/ncov/open/africa/6m?f_region=Africa)                        | [africa/all-time](/ncov/open/africa/all-time?f_region=Africa)                        |
| **Asia**          | [asia/2m](/ncov/open/asia/2m?f_region=Asia)                              | [asia/6m](/ncov/open/asia/6m?f_region=Asia)                              | [asia/all-time](/ncov/open/asia/all-time?f_region=Asia)                              |
| **Europe**        | [europe/2m](/ncov/open/europe/2m?f_region=Europe)                        | [europe/6m](/ncov/open/europe/6m?f_region=Europe)                        | [europe/all-time](/ncov/open/europe/all-time?f_region=Europe)                        |
| **North America** | [north-america/2m](/ncov/open/north-america/2m?f_region=North%20America) | [north-america/6m](/ncov/open/north-america/6m?f_region=North%20America) | [north-america/all-time](/ncov/open/north-america/all-time?f_region=North%20America) |
| **Oceania**       | [oceania/2m](/ncov/open/oceania/2m?f_region=Oceania)                     | [oceania/6m](/ncov/open/oceania/6m?f_region=Oceania)                     | [oceania/all-time](/ncov/open/oceania/all-time?f_region=Oceania)                     |
| **South America** | [south-america/2m](/ncov/open/south-america/2m?f_region=South%20America) | [south-america/6m](/ncov/open/south-america/6m?f_region=South%20America) | [south-america/all-time](/ncov/open/south-america/all-time?f_region=South%20America) |

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of 8 &times; 10^-4 subs per site per year. Mutational fitness is calculated using results from [Obermeyer et al (under review)](https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1). Full details on bioinformatic processing can be found [here](https://github.com/nextstrain/ncov).

The analysis on this page uses data from NCBI GenBank as a source following [Open Data principles](https://opendatahandbook.org/guide/en/what-is-open-data/), such that we can make input data and intermediate files available for further analysis. Open Data is data that can be freely used, re-used and redistributed by anyone - subject only, at most, to the requirement to attribute and sharealike. But be aware that not all regions are well represented in open databases and some of the above trees might lack recent data from particular geographic regions.

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work in open databases. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain. An attribution table is available by clicking on "Download Data" at the bottom of the page and then clicking on "Strain Metadata" in the resulting dialog box.

To maximize the utility and visibility of these generously shared data, [we provide preprocessed files that can serve as a starting point for additional analyses](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html).

#### All sequences and metadata

- [metadata.tsv.gz](https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz)
- [sequences.fasta.xz](https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz)
- [aligned.fasta.xz](https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz)

Now also available with `zstd` compression allowing much faster decompression:

- [metadata.tsv.zst](https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst)
- [sequences.fasta.zst](https://data.nextstrain.org/files/ncov/open/sequences.fasta.zst)
- [aligned.fasta.zst](https://data.nextstrain.org/files/ncov/open/aligned.fasta.zst)

#### Subsampled sequences and intermediate files

The files below exist for every region (`global`, `africa`, `asia`, `europe`, `north-america`, `oceania` and `south-america`) and correspond to each region's 6 month timespan build (e.g. `global/6m`, `africa/6m`, `asia/6m`, etc).
Files for the `2m` and `all-time` builds (e.g. `global/2m`, `global/all-time`, etc.) are not yet available.
The links below refer to the `${BUILD_PART_0}` region; substitute `${BUILD_PART_0}` with another region name in the links if desired.

- [${BUILD_PART_0}/6m metadata.tsv.xz](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/metadata.tsv.xz)
- [${BUILD_PART_0}/6m sequences.fasta.xz](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/sequences.fasta.xz)
- [${BUILD_PART_0}/6m aligned.fasta.xz](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/aligned.fasta.xz)
- [${BUILD_PART_0}/6m Auspice tree](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/${BUILD_PART_0}.json)
- [${BUILD_PART_0}/6m Auspice root sequence](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/${BUILD_PART_0}_root-sequence.json)
- [${BUILD_PART_0}/6m Auspice tip frequencies](https://data.nextstrain.org/files/ncov/open/${BUILD_PART_0}/${BUILD_PART_0}_tip-frequencies.json)
