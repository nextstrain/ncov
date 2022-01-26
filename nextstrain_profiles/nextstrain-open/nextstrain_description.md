Compiled Nextstrain SARS-CoV-2 resources are available at [nextstrain.org/sars-cov-2](https://nextstrain.org/sars-cov-2/). Follow [@nextstrain](https://twitter.com/nextstrain) for continual data updates.

This phylogeny shows evolutionary relationships of SARS-CoV-2 viruses from the ongoing COVID-19 pandemic. Although the genetic relationships among sampled viruses are quite clear, there is considerable uncertainty surrounding estimates of specific transmission dates and in reconstruction of geographic spread. Please be aware that specific inferred geographic transmission patterns and temporal estimates are only a hypothesis.

There are hundreds of thousands of complete SARS-CoV-2 genomes available on open databases and this number increases every day, but geographical representation varies. This visualization can only handle ~3000 genomes in a single view for performance and legibility reasons. Because of this we subsample available genome data for these analysis views. Our primary [global analysis](/ncov/open/global/) subsamples to ~600 genomes per continental region with ~400 from the previous 4 months and ~200 from before this. This results in a more equitable global sequence distribution, but hides samples available from regions that are doing lots of sequencing. To mitigate against this, we've set up separate analyses to focus on particular regions. They are available on the "Dataset" dropdown on the left or by clicking on the following links: [Africa](/ncov/open/africa?f_region=Africa), [Asia](/ncov/open/asia?f_region=Asia), [Europe](/ncov/open/europe?f_region=Europe), [North America](/ncov/open/north-america?f_region=North%20America), [Oceania](/ncov/open/oceania?f_region=Oceania) and [South America](/ncov/open/south-america?f_region=South%20America).

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of 8 &times; 10^-4 subs per site per year. Mutational fitness is calculated using results from [Obermeyer et al (under review)](https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1). Full details on bioinformatic processing can be found [here](https://github.com/nextstrain/ncov).

The analysis on this page uses data from NCBI GenBank as a source following Open Data principles, such that we can make input data and intermediate files available for further analysis (see below). Open Data is data that can be freely used, re-used and redistributed by anyone - subject only, at most, to the requirement to attribute and sharealike. But be aware that not all regions are well represented in open databases and some of the above trees might lack recent data from particular geographic regions.

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work in open databases. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain. An attribution table is available by clicking on "Download Data" at the bottom of the page and then clicking on "Strain Metadata" in the resulting dialog box.

To maximize the utility and visibility of these generously shared data, [we provide preprocessed files that can serve as a starting point for additional analyses](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html).

### All sequences and metadata

#### Ingested and parsed data

 * [sequences.fasta.xz](https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz)
 * [metadata.tsv.gz](https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz)

#### Pre-processed files

 * [aligned.fasta.xz](https://data.nextstrain.org/files/ncov/open/aligned.fasta.xz)
 * [filtered.fasta.xz](https://data.nextstrain.org/files/ncov/open/filtered.fasta.xz)
 * [mutation-summary.tsv.xz](https://data.nextstrain.org/files/ncov/open/mutation-summary.tsv.xz)

### Subsampled sequences and intermediate files

The files below exist for the `global` and the regional builds (`africa`, `asia`, `europe`, `north-america`, `oceania` and `south-america`).
The links below refer to the `${BUILD}` build, substitute `${BUILD}` with another build name in the links if desired.

 * [${BUILD}/sequences.fasta.xz](https://data.nextstrain.org/files/ncov/open/${BUILD}/sequences.fasta.xz)
 * [${BUILD}/metadata.tsv.xz](https://data.nextstrain.org/files/ncov/open/${BUILD}/metadata.tsv.xz)
 * [${BUILD}/aligned.fasta.xz](https://data.nextstrain.org/files/ncov/open/${BUILD}/aligned.fasta.xz)
 * [${BUILD} auspice tree](https://data.nextstrain.org/files/ncov/open/${BUILD}/${BUILD}.json)
 * [${BUILD} auspice root sequence](https://data.nextstrain.org/files/ncov/open/${BUILD}/${BUILD}_root-sequence.json)
 * [${BUILD} auspice tip frequencies](https://data.nextstrain.org/files/ncov/open/${BUILD}/${BUILD}_tip-frequencies.json)
