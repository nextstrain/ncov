# ncov_wa

## Build Overview
- **Build Name**: [ncov_wa_six_mon]
- **Pathogen/Strain**: [ncov]
- **Scope**: [WGS of SARS-CoV-2 in Washington state]
- **Purpose**: [Genomic surveillance of SARS-CoV-2 in Washington State for past six months]
- **Nextstrain Build Location**: [Washington-focused SARS-CoV-2 genomic analysis: Past six months](https://nextstrain.org/groups/waphl/ncov/wa/6m)

## Table of Contents:
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

## Getting Started

### Data Sources & Inputs
This build uses the Full Remote Dataset and Global Remote Datasets available on [Nextstrain](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html). This build is designed to pull Washington state sequences and metadata from the full remote dataset as the inputs to the ncov Nextstrain pipeline. The Global dataset (alignment and metadata) is used for contextual sequences in the build. To include more contextualization, one could use the Full Remote Dataset for the contextual sequences, however doing so may require AWS Batch to subsample from the dataset.

### Setup & dependencies
See installation.

### Installation
First, install the [ncov nextstrain pipeline](https://github.com/nextstrain/ncov) and clone the ncov repository using `git clone https://github.com/nextstrain/ncov` or `gh repo clone nextstrain/ncov`.

### Clone the Repository
Clone this repository in the `ncov` folder. You can do this in the command-line terminal by navigating to the `ncov` repository using `cd ncov` and then cloning the repository using `git clone https://github.com/DOH-SML1303/ncov_wa.git` or `gh repo clone DOH-SML1303/ncov_wa`.

## Run the Build
This ncov Nexstrain build sources data from Genbank and includes a 6m build. If you're running Nextstrain in a conda environment or `Nextstrain shell` then you want to make sure you pull the latest ncov github repository updates first by running `git pull` in the `ncov` directory, activating the conda environment using `conda activate nextstrain` or Nexstrain shell using `Nextstrain shell .` followed by `nextstrain update` to update Nextstrain. (To update the `Nextstrain shell`, you must run `nextstrain update` outside of the shell) It's recommended to pull updates prior to running the pipeline. The same could also be said for this repo as well! :)

### To run the builds using inputs stored on an AWS Bucket:
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
`nextstrain build --aws-batch-s3-bucket bucket-name --cpus=6 . --configfile ncov_wa/config/builds.yaml`

### Run the builds locally
`nextstrain build --cpus=6 . --configfile ncov_wa/config/builds.yaml`

## Visualizing the results
You can check your results once the pipeline is done running using `nextstrain view auspice`

## Repository FIle Structure Overview
The file hierachy for this customized build:
```
ncov_wa/
├──config/
|   ├──auspice_config.json         #variables to include in Color By feature
|   ├──builds.yaml                 #the builds file that customizes nextstrain build
|   ├──colors.tsv                  #WA county colors
|   ├──config.yaml                 #file that includes dependencies and path to the builds.yaml
|   ├──description.md
├──data/
|   ├──county_metadata.tsv         #to add WA counties to the metadata so they can be included in the build
|   ├──headers.tsv                 #needed for the smk workflow to create metadata file
├──scripts/
|   ├──wa-nextstrain-update-location-genbank.py #adds county metadata to the filtered wa seqs metadata
|   ├──filter_wa_metadata.sh       #for the smk workflow to pull the WA metadata from the full remote dataset
|   ├──filter_wa_sequences.sh      #for the smk workflow to pull the WA sequences from the full remote dataset
|   ├──pull_full_data.sh           #for the smk workflow to pull the full remote dataset to filter out anything that's not WA seqs and metadata
├──workflow/
|   ├──filter_wa_data.smk          #pulls the full data and then filters for WA data to be the input into the Nextstrain build
```

## Files that may need to be changed
When you pull updates for the ncov repo there are a few files that you want to keep an eye for for any changes. This includes the following files the default ncov build:
- `ncov/defaults/auspice_config.json`
- `ncov/nextstrain_profiles/.../builds.yaml`

If there are any changes to these two files then changes may need to be made to their custom counterparts in this focused build.
- Changes to `ncov/defaults/auspice_config.json` > make changes to > `ncov_wa/config/auspice_config.json`
- Changes to `ncov/nextstrain_profiles/.../builds.yaml` > *may require changes to* > `ncov_wa/config/builds.yaml`

## Expected Outputs
- `ncov/auspice/ncov_ncov_wa_six_mon.json`
- `ncov/auspice/ncov_ncov_wa_six_mon_root-sequence.json`
- `ncov/auspice/ncov_ncov_wa_six_mon_tip-frequencies.json`

## Scientific Decisions

## Adapting for Another Jursidiction

## Contributing

## License

## Acknowledgments
