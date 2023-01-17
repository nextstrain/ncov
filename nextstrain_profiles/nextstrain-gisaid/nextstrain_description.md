Compiled Nextstrain SARS-CoV-2 resources are available at [nextstrain.org/sars-cov-2](https://nextstrain.org/sars-cov-2/). Follow [@nextstrain](https://twitter.com/nextstrain) for updates.

This phylogeny shows evolutionary relationships of SARS-CoV-2 viruses from the ongoing COVID-19 pandemic. Although the genetic relationships among sampled viruses are generally quite clear, there is considerable uncertainty surrounding estimates of specific transmission dates and in reconstruction of geographic spread. Please be aware that specific inferred geographic transmission patterns and temporal estimates are only a hypothesis.

There are millions of complete SARS-CoV-2 genomes available and this number increases every day. This visualization can only handle ~4000 genomes in a single view for performance and legibility reasons. Because of this we subsample available genome data for our analysis views. We provision multiple views to focus subsampling on different geographic regions and different time periods. These views are available through the "Dataset" dropdown on the left or by clicking on the following links:

&nbsp;            | past 1 month                                                               | past 2 months                                                              | past 6 months                                                              | all time
----------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------
**global**        | [global/1m](/ncov/gisaid/global/1m)                                        | [global/2m](/ncov/gisaid/global/2m)                                        | [global/6m](/ncov/gisaid/global/6m)                                        | [global/all-time](/ncov/gisaid/global/all-time)                                        |
**Africa**        | [africa/1m](/ncov/gisaid/africa/1m?f_region=Africa)                        | [africa/2m](/ncov/gisaid/africa/2m?f_region=Africa)                        | [africa/6m](/ncov/gisaid/africa/6m?f_region=Africa)                        | [africa/all-time](/ncov/gisaid/africa/all-time?f_region=Africa)                        |
**Asia**          | [asia/1m](/ncov/gisaid/asia/1m?f_region=Asia)                              | [asia/2m](/ncov/gisaid/asia/2m?f_region=Asia)                              | [asia/6m](/ncov/gisaid/asia/6m?f_region=Asia)                              | [asia/all-time](/ncov/gisaid/asia/all-time?f_region=Asia)                              |
**Europe**        | [europe/1m](/ncov/gisaid/europe/1m?f_region=Europe)                        | [europe/2m](/ncov/gisaid/europe/2m?f_region=Europe)                        | [europe/6m](/ncov/gisaid/europe/6m?f_region=Europe)                        | [europe/all-time](/ncov/gisaid/europe/all-time?f_region=Europe)                        |
**North America** | [north-america/1m](/ncov/gisaid/north-america/1m?f_region=North%20America) | [north-america/2m](/ncov/gisaid/north-america/2m?f_region=North%20America) | [north-america/6m](/ncov/gisaid/north-america/6m?f_region=North%20America) | [north-america/all-time](/ncov/gisaid/north-america/all-time?f_region=North%20America) |
**Oceania**       | [oceania/1m](/ncov/gisaid/oceania/1m?f_region=Oceania)                     | [oceania/2m](/ncov/gisaid/oceania/2m?f_region=Oceania)                     | [oceania/6m](/ncov/gisaid/oceania/6m?f_region=Oceania)                     | [oceania/all-time](/ncov/gisaid/oceania/all-time?f_region=Oceania)                     |
**South America** | [south-america/1m](/ncov/gisaid/south-america/1m?f_region=South%20America) | [south-america/2m](/ncov/gisaid/south-america/2m?f_region=South%20America) | [south-america/6m](/ncov/gisaid/south-america/6m?f_region=South%20America) | [south-america/all-time](/ncov/gisaid/south-america/all-time?f_region=South%20America) |

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of 8 &times; 10^-4 subs per site per year. Mutational fitness is calculated using results from [Obermeyer et al (under review)](https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1). Full details on bioinformatic processing can be found [here](https://github.com/nextstrain/ncov).

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata made available through [GISAID](https://gisaid.org) on which this research is based. An attribution table is available by clicking on "Download Data" at the bottom of the page and then clicking on "Acknowledgments" in the resulting dialog box.

At the specific request of GISAID, we:
 - maintain the prefix `hCoV-19/` in the names of viral isolates
 - disable download of full metadata TSV and provide instead an acknowledgments TSV in the "Download Data" link at the bottom of the page
 - refrain from sharing alignments or other intermediate files computed in our pipeline

#### Reusing code or images

All source code for Auspice, the visualization tool, is freely available under the terms of the [GNU Affero General Public License](https://github.com/nextstrain/auspice/blob/HEAD/LICENSE.txt).

Screenshots may be used under a [CC-BY-4.0 license](https://creativecommons.org/licenses/by/4.0/) and attribution to nextstrain.org must be provided. A high-quality download option is available by clicking the **DOWNLOAD DATA** button at the bottom of the page and selecting **SCREENSHOT (SVG)**.
