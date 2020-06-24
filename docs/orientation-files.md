
## Overview of this repository (i.e., what do these files do?)

The files in this repository fall into one of these categories:  
* Input files  
* Output files and directories  
* Workflow configuration files we might want to customize  
* Workflow configuration files we don't need to touch  
* Documentation  

We'll walk through all of the files one by one, but here are the most important ones for your reference:  

|Category| Directory | File | Description | Configuration|  
|-----|-----|-----|------|-----|
|Input file|`./input_data/`|`sequences.fasta`|**Genomic sequences; IDs must match `strain` column in `metadata.tsv`**| See ['Preparing your data'](data-prep.md)
|Input file|`./input_data/`|`metadata.tsv`|**Tab-delimited description of strain (i.e., sample) attributes**|See ['Preparing your data'](data-prep.md)|
|Output file|`./auspice/`|`buildName.json`|**Output file for visualization in auspice**||
|Customizable workflow file|`./my_config/<mybuildname>/`|`builds.yaml`|**Specify and define all the builds you'd like to run**|See our [customization guide](customizing-analysis.md)|
|Customizable workflow file|`./my_config/<mybuildname>/`|`config.yaml`|**Analysis configuration file; parameterize your analyses here**|See our [customization guide](customizing-analysis.md)|


-----


## Input files  

| Directory | File | Description | Configuration|  
|-----|-----|-----|------|
|`./input_data/`|`sequences.fasta`|**Genomic sequences; IDs must match `strain` column in `metadata.tsv`**| See ['Preparing your data'](data-prep.md)
|`./input_data/`|`metadata.tsv`|**Tab-delimited description of strain (i.e., sample) attributes**|See ['Preparing your data'](data-prep.md)|
|`./default_config/`|`include.txt`| List of strain names to _include_ during subsampling and filtering | One strain name per line|  
|`./default_config/`|`exclude.txt`|List of strain names to _exclude_ during subsampling and filtering|One strain name per line|


## Output files and directories  

| Directory | File | Description |
|-----|-----|-----|
|`./auspice/`|`buildName.json`|**Output file for visualization in auspice**|
|`./intermediate_results/`|`aligned.fasta`, `sequence-disagnostics.tsv`, etc.|Raw results files (dependencies) that are shared across all `builds`|
|`./intermediate_results/<buildName>/`|`tree.nwk`, `aa_mutations.json`, etc.|Raw results files (dependencies) that are specific to a single `build`|
|`./logs/`|`.log` files|Error messages and other information about the run|


## Workflow configuration files we might want to customize  

| Directory | File | Description | Configuration |
|-----|-----|-----|----|
|`./my_config/<mybuildname>/`|`builds.yaml`|**Specify and define all the builds you'd like to run**|See our [customization guide](customizing-analysis.md)|
|`./my_config/<mybuildname>/`|`config.yaml`|**Analysis configuration file; parameterize your analyses here**|See our [customization guide](customizing-analysis.md)|
|`./default_config/`|`default_config.yaml`|**Default analysis configuration file**|Override in `./my_config/.../config.yaml`|
|`./default_config/`|`default_auspice_config.json`|**Default visualization configuration file**|Override in `./my_config/.../auspice_config.yaml`|


## Workflow configuration files we don't need to touch  

| Directory | File | Description | Configuration|
|-----|-----|-----|-----|
|`./`|`Snakefile`|Entry point for `snakemake` commands; validates input.|No modification needed|
|`./workflow/snakemake_rules/`|`main_workflow.smk`|Defines rules for running each step in the analysis|Modify your `config` file, rather than hardcode changes into the snakemake file itself|
|`./workflow/envs/`|`nextstrain.yaml`|Specifies computing environment needed to run workflow with the `--use-conda` flag|No modification needed|
|`./workflow/schemas/`|`config.schema.yaml`|Defines format (e.g., required fields and types) for  `config.yaml` files.|Useful reference, but no modification needed.|
|`./scripts/`| `add_priorities_to_meta.py`, etc.| Helper scripts for common tasks | No modification needed |

## [Previous Section: Orientation: analysis workflow](orientation-workflow.md)
## [Next Section: Running & troubleshooting](running.md)
