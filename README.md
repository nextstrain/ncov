# About  

This repository analyzes viral genomes using [Nextstrain](https://nextstrain.org) to understand how SARS-CoV-2, the virus that is responsible for the COVID-19 pandemic, evolves and spreads.

We maintain a number of publicly-available builds, visible at [nextstrain.org/ncov](https://nextstrain.org/ncov).

---
# Resources

### Use Nextstrain to analyze your SARS-CoV-2 data  

**We've written a comprehensive guide to get you up and running in <1 hr. Click on the below links to follow it. It covers:**
* [**Preparing your data**](./docs/data-prep.md) _(Start here)_
* [**Setup and installation**](./docs/setup.md)
* [**Orientation: analysis workflow**](.docs/orientation-workflow.md)
* [**Orientation: which files should I touch?**](.docs/orientation-files.md)
* [**Running & troubleshooting**](.docs/running.md)
* [**Customizing your analysis**](.docs/customizing-analysis.md)
* [**Customizing your visualization**](.docs/customizing-visualization.md)
* [**Options for visualizing and sharing results**](.docs/sharing.md) (including working with sensitive metadata)  
* [**Interpreting your results**](.docs/interpretation.md)
* [**Writing a narrative to highlight key findings**](.docs/narratives.md)


### Download formatted datasets  

The hCoV-19 / SARS-CoV-2 genomes were generously shared via GISAID. We gratefully acknowledge the Authors, Originating and Submitting laboratories of the genetic sequence and metadata made available through GISAID on which this research is based.

In order to download the GISAID data to run the analysis yourself, please see [this guide](XXX).  
> Please note that `input_data/metadata.tsv` is no longer included as part of this repo. However, we provide continually-updated, pre-formatted metadata & fasta files for download through GISAID.

### Read previous Situation Reports  
We issued weekly Situation Reports for the first ~5 months of the pandemic. You can find the Reports and their translations [here](https://nextstrain.org/ncov-sit-reps).

### FAQs  

- Can't find your sequences in Nextstrain? Check [here](./docs/data_faq.md) for common reasons why your sequences may not be appearing.
- For information about how clades are defined, and the currently named clades, please see [here](./contributor_docs/clades.md).

### Bioinformatics notes

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of [8 &times; 10^-4 subs per site per year](http://virological.org/t/phylodynamic-analysis-176-genomes-6-mar-2020/356). There were SNPs present in the nCoV samples in the first and last few bases of the alignment that were masked as likely sequencing artifacts.

---

# Contributing

We welcome contributions from the community! Please note that we strictly adhere to the [Contributor Covenant Code of Conduct](https://github.com/nextstrain/.github/blob/master/CODE_OF_CONDUCT.md).

### Contributing to software or documentation   
Please see our [Contributor Guide](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md) to get started!

### Contributing data  
**Please note that we automatically pick up any SARS-CoV-2 data that is submitted to GISAID.**  

If you're a lab and you'd like to get started sequencing, please see:  
* [Protocols from the ARTIC network](https://www.protocols.io/groups/artic/publications)  
* [Funding opportunities for sequencing efforts](https://twitter.com/firefoxx66/status/1242147905768751106)  
* Or, if these don't meet your needs, [get in touch](mailto:hello@nextstrain.org)  

---

# Get in touch  

To report a bug, error, or feature request, please [open an isssue](https://github.com/nextstrain/ncov/issues).  

For questions, head over to the [discussion board](https://discussion.nextstrain.org/); we're happy to help!
