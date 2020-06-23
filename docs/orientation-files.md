
## Overview of this repository (i.e., what do these files do?)

The files in this repository fall into one of XXX categories:  
* Input files  
* Output files and directories  
* Workflow configuration files we might want to customize  
* Workflow configuration files we don't need to touch  
* Documentation  

We'll walk through these one by one; **the 5 most important files are bolded**.

## Input files  

<!-- XXX make file names into links -->
| Directory | File | Description | Configuration|  
|-----|-----|-----|------|
|`./data/`|`sequences.fasta`|**Genomic sequences; IDs must match `strain` column in `metadata.tsv`**| See ['Preparing your data'](XXX)
|`./data/`|`metadata.tsv`|**Tab-delimited description of strain (i.e., sample) attributes**|See ['Preparing your data'](XXX)|
|`./config/`|`include.txt`| List of strain names to forcibly _include_ during subsampling and filtering | One strain name per line|  
|`./config/`|`exclude.txt`|List of strain names to forcibly _exclude_ during subsampling and filtering|One strain name per line|


## Output files and directories  

| Directory | File | Description |
|-----|-----|-----|
|`./auspice/`|`buildName.json`|**Output file for visualization in auspice**|
|`./results/`|`aligned.fasta`, `sequence-disagnostics.tsv`, etc.|Raw results files (dependencies) that are shared across all `builds`|
|`./results/<buildName>/`|`tree.nwk`, `aa_mutations.json`, etc.|Raw results files (dependencies) that are specific to a single `build`|
|`./logs/`|`.log` files|Error messages and other information about the run|


## Workflow configuration files we might want to customize  

| Directory | File | Description | Configuration |
|-----|-----|-----|----|
|`./my_config/<mybuildname>/mybuilds.yaml`|**Specify and define all the builds you'd like to run**|See our [customization guide](XXX)|
|`./my_config/<mybuildname>/myconfig.yaml`|**Analysis configuration file; parameterize your analyses here**|See our [customization guide](XXX)|
|`./default_config/`|`default_config.yaml`|**Default analysis configuration file**|Override these settings in `./my_config/.../config.yaml`|
|`./default_config/`|`default_auspice_config.json`|**Default visualization configuration file**|Override these settings in `./my_config/.../auspice_config.yaml`|XXX|


## Workflow configuration files we don't need to touch  
| Directory | File | Description | Configuration|
|-----|-----|-----|-----|
|`./`|`Snakefile`|Entry point for `snakemake` commands; validates input.|No modification needed|
|`./rules/`|`builds.smk`|Defines rules for running each step in the analysis|Modify your `config` file, rather than hardcode changes into the snakemake file itself|
|`./envs/`|`nextstrain.yaml`|Specifies computing environment needed to run workflow with the `--use-conda` flag|No modification needed|
|`./schemas/`|`config.schema.yaml`|Defines format (e.g., required fields and types) for  `config.yaml` files.|Useful reference, but no modification needed.|
|`./scripts/`| `add_priorities_to_meta.py`, etc.| Helper scripts for common tasks | No modification needed |

## [Previous Section: Orientation: analysis workflow](./docs/orientation-workflow.md)
## [Next Section: Orientation: Running & troubleshooting](./docs/running.md)