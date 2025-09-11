## SARS-CoV-2 Washington focused build

### Build Overview
- **Build Name**: SARS-CoV-2 Washington focused build
- **Pathogen/Strain**: SARS-CoV-2
- **Scope**: Whole Genome Sequences of SARS-CoV-2 in Washington state from the past year
- **Purpose**: This repository contains the Nextstrain build for the genomic surveillance of SARS-CoV-2 in Washington State for past year.
- **Nextstrain Build Location**: [Washington-focused SARS-CoV-2 genomic analysis: Past year](https://nextstrain.org/groups/waphl/ncov/wa/1y)

### Table of Contents:
- [Getting Started](#getting-started)
  - [Data Sources & Inputs](#data-sources--inputs)
  - [Setup & Dependencies](#setup--dependencies)
    - [Installation](#installation)
    - [Clone the repository](#clone-the-repository)
- [Run the Build](#run-the-build)
  - [To run the builds using inputs stored on an AWS Bucket](#to-run-the-builds-using-inputs-stored-on-an-AWS-Bucket)
  - [Run the builds locally](#run-the-builds-locally)
- [Visualizing the results](#Visualizing-the-resules)
- [Repository File Structure Overview](#repository-file-structure-overview)
- [Expected Outputs](#expected-outputs)
- [Scientific Decisions](#scientific-decisions)
- [Adapting for Another Jurisdiction](#adapting-for-another-jurisdiction)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

### Pathogen Epidemiology

Overview:

- Taxonomic designations
- Geographic distribution and seasonality
- Public Health importance
- Genomic relevance
- Additional resources

### Getting Started
This build utilizes the [Nextstrain.org remote datasets](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html) to produce a Washington-focused SC2 Nextstrain build that can be used for genomic surveillance purposes.

Some high-level build features and capabilities are:
- **1 year Washington focus sampling**: All Washington sequences from the last year are included in this build.
- **Tiered subsampling**: Additional sequences from the rest of the USA & the world are selected by genetic similarity to the state-level sequences. Additionally, earlier sequences from Washington and globally are provided for temporal context.

### Data Sources & Inputs
This build uses NCBI data and the SARS-Cov-2 Global Remote Dataset available on [Nextstrain.org](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html). The Remote Dataset data is sourced from GenBank cleaned/maintainted by the Nextstrain team.  This build pulls in subsets Washington State sequences  and metadata from GenBank, and pulls in contextual data from Nextstrain Global Remote Dataset that are the inputs to the ncov Nextstrain pipeline.

To include more contextualization, one could use the Full SARS-Covo2 Remote Dataset for the contextual sequences, however doing so may require AWS Batch to subsample from the dataset.

- **Sequence Data**: GenBank SARS-Cov-2 data from Datasets and Nextstrain.org SC2 Remote Dataset sourced GenBank
- **Metadata**: GenBank SARS-Cov-2 data from Datasets, Nextstrain.org SC2 Remote Dataset sourced GenBank and WA DOH county-level data
- **Expected Inputs**:
    - `ncov_wa/data/county_metadata.csv` (contains most recent line list of GenBank accession number and Washington State county designation)
    -  Other sequencing and metadata will be automatically downloaded and ingested as part of the pipeline

### Setup & dependencies
#### Installation
Ensure that you have [Nextstrain](https://docs.nextstrain.org/en/latest/install.html) installed.

To check that Nextstrain is installed:
```
nextstrain check-setup
```
If Nextstrain is not installed, follow [Nextstrain installation guidelines](https://docs.nextstrain.org/en/latest/install.html)

#### Clone this ncov repository:
Clone this repository by running:

```
git clone https://github.com/NW-PaGe/ncov.git
```


<!--### To run the builds using inputs stored on an AWS Bucket:
You can configure your `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` in your AWS credentials file which can be accessed in terminal using `nano ~/.aws/credentials`, or you can simply export the environmental variables upon opening a terminal window using:
`export AWS_ACCESS_KEY_ID=
export AWS_SECRET_ACCESS_KEY=`

There's some additional modifications you would have to include in the `ncov_wa/config/builds.yaml` to ensure the pipeline know to read from your bucket. You could just include the following code at the top of the file:

```
S3_DST_BUCKET: <bucket path>
S3_DST_COMPRESSION: "xz" #if your outputs are compressed
S3_DST_ORIGINS: [ncov-wa] #name of your inputs
upload:
  - build-files
```

If you're running Batch then you need to make sure all of the information is included in your `~/.nextstrain/config`. File. See [this documentation](https://docs.nextstrain.org/projects/cli/en/stable/aws-batch/) for more information.

To run the builds with your data stored in an AWS Bucket, navigate to the `ncov` directory and run:
`nextstrain build --aws-batch-s3-bucket bucket-name --cpus=6 . --configfile ncov_wa/config/builds.yaml` -->

### Run the builds locally

To run the build, make sure you are in the correct directory file "ncov".  The below code specifies how many CPUs to use as well as which config file to use. In this case, we are specifying to use the  **ncov_wa/config/build.yaml** with our Washington-specific parameters.

```
nextstrain build --cpus=6 . --configfile ncov_wa/config/builds.yaml
```

### Repository File Structure Overview
The file structure of the repository is as follows with `*`  denoting folders that are the build's expected outputs.

```
.
├──ncov_wa/
|   ├──config/
|   |   ├──auspice_config.json         #  Variables to include in Color-By feature
|   |   ├──builds.yaml                 #  The builds file that customizes the Nextstrain build; includes subsampling scheme, filters and sites to be masked
|   |   ├──colors.tsv                  #  File the designates the colors for division, region, and location in Auspice
|   |   ├──config.yaml                 #  File that includes dependencies and path to the builds.yaml
|   |   ├──description.md
|   ├──data/
|   |   ├──county_metadata.csv         #  Most recent CSV line list of  GenBank accession and associated Washington County. This file need to be updated before running build.
|   |   ├──headers.tsv                 #  TSV file with the headers of interest for the metadata file.  This file is needed for the smk workflow to use the metadata file for Washington filtering.
|   |   ├──tmp.tsv                     # Temporary file used in filtering Washington sequences step
|   ├──scripts/
|   |   ├──filter_wa_metadata.sh                            # Bash script that filters the remote metadata set to WA metadata. This is called in the workflow/filter_wa_datas.smk file
|   |   ├──filter_wa_sequences.sh                           # Bash script that filters the remote dataset sequences to WA sequences. This is called in the workflow/filter_wa_data.smk file
|   |   ├──pull_full_data.sh                                # Bash script that pulls in the full remote dataset. This is called in the workflow/filter_wa_data.smk file
|   |   ├──wa-nextstrain-update-location-genbank.py         # Python script that joins the county_metadata.csv file to the WA metadata
|   ├──workflow/
|   |   ├──filter_wa_data.smk                              # Pulls the full remote dataset and then filters for WA data to be the input into the Nextstrain build
|   |   ├──add_to_builds.smk                               # Calculates dates and specifies date scheme for timeframe of data the build should include
```
More details on the file structure of this build can be found [here](https://github.com/NW-PaGe/ncov/wiki/_new)

### Visualizing the results
You can check your results once the pipeline is done running using `nextstrain view auspice` or by uploading the output JSON file into [auspice.us](https://auspice.us/).

### Files that  need to be updated
When running the build, the *county_metadata.csv* should be updated to capture the most up-to-date county data. This metadata file is generated by WA DOH and contains two columns: **SEQUENCE_GENBANK_STRAIN** containing GenkBank accession IDs that match to the sequence FASTA headers, and **COUNTY_NAME** column listing the associated county for each sequence.

<!-- When you pull updates for the ncov repo there are a few files that you want to keep an eye for for any changes. This includes the following files the default ncov build:
- `ncov/defaults/auspice_config.json`
- `ncov/nextstrain_profiles/.../builds.yaml` <--

If there are any changes to these two files then changes may need to be made to their custom counterparts in this focused build.
- Changes to `ncov/defaults/auspice_config.json` > make changes to > `ncov_wa/config/auspice_config.json`
- Changes to `ncov/nextstrain_profiles/.../builds.yaml` > *may require changes to* > `ncov_wa/config/builds.yaml` -->

### Expected Outputs

Within the `ncov/auspice/` folder, the expected outputs include:
- `ncov_ncov_wa_1y.json`
- `ncov_ncov_wa_1y_root-sequence.json`
- `ncov_ncov_wa_1y_tip-frequencies.json`

### Scientific Decisions
- **Subsampling**:
  - **1 year Washington focus sampling**: Subsampling includes all Washington sequences (no maximum number of sequences) from the past year
  - **Contextual proximity sampling**: Subsampling includes 1000 sequences sampled from 2020 through current. This sampling helps to accurately reconstruct the number of introduction.  Proximity sampling selects sequences as close as possible to the focal samples (Currently set to Washington).  The genetic proximity between sequences in the focal set to other sequences are calculated in the [priorities.py](https://github.com/nextstrain/ncov/blob/5555ece97bafe1aa2cb19dcaac183d5a718d29fa/scripts/priorities.py) script.
    - **Crowding penalty**: The crowding penalty in proximity subsampling controls how strongly the subsampling penalizes sequences that are genetically similar to each other. This build set the crowding penalty to 0. The default setting is 0.25.  A crowding penalty value closer to 1 creates a bushier tree and discourages sequence redundancy. A crowding penalty closer to 0 allows more clustering. A crowding penalty of 0 disables crowding.
  - **Contextual random sampling**: Subsampling includes 500 sequences sampled over month-year that allow for accurate clade timing in the tree.
- **Reference selection**: [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) is used as the reference because it is the complete genome of the SARS-CoV-2 Wuhan strain collected in December 2019.
- **Clade labeling**: Internal clade labels are included in the tree through the [main_workflow.smk](https://github.com/nextstrain/ncov/blob/5555ece97bafe1aa2cb19dcaac183d5a718d29fa/workflow/snakemake_rules/main_workflow.smk#L954)

## Adapting for Another Jursidiction
- The jurisdiction-focused sampling time frame of the build can be changed. It is currently set up to focus on the last year of Washington sequences, but this time frame can be altered to be shorter/longer by adjusting the add_to_builds.smk and the build.yaml subsampling scheme.
- To adapt the build to a new jurisdiction, the current filters for Washington should be changed to filter for jurisdiction of interest.  These filtering steps are in the [filter_wa_metadata.sh](https://github.com/NW-PaGe/ncov_wa/blob/main/scripts/filter_wa_metadata.sh) bash script that pattern matches the metadata, and that is called within the [filter_wa.smk](https://github.com/NW-PaGe/ncov_wa/blob/main/workflow/filter_wa_data.smk) workflow. Note: when working with bash scripts, be careful about editing the files in a Windows application, and be sure the files are saved with only the line feed character (LF) instead of the carriage return plus line feed (CRLF).
- The [colors.tsv](https://github.com/NW-PaGe/ncov_wa/blob/main/config/colors.tsv) file can be adapted to change colors visualized in Auspice Color-By. The tsv file should include the divisions of interest that are to appear in the Color-By.

## Contributing
For any questions please submit them to our [Discussions](https://github.com/NW-PaGe/ncov_wa/discussions) page otherwise software issues and requests can be logged as a Git [Issue](https://github.com/NW-PaGe/ncov_wa/issues).
## License
This project is licensed under a modified GPL-3.0 License.
You may use, modify, and distribute this work, but commercial use is strictly prohibited without prior written permission.

## License

## Acknowledgments
