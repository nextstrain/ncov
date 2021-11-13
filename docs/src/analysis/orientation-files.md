# Overview of this repository (i.e., what do these files do?)

The files in this repository fall into one of these categories:
* Input files
* Output files and directories
* Workflow configuration files we might want to customize
* Workflow configuration files we don't need to touch
* Documentation

We'll walk through all of the files one by one, but here are the most important ones for your reference:

  * Input files
    * `data/metadata.tsv`: tab-delimited description of strain (i.e., sample) attributes
    * `data/sequences.fasta`: genomic sequences whose ids must match the `strain` column in `metadata.tsv`. [See the data preparation guide](data-prep.md).
    * `my_profiles/<your_profile>/builds.yaml`: workflow configuration file where you can define and parameterize the builds you want to run. The directory name `your_profile` is the name of your custom analysis profile where you store this configuration and other custom files for the analysis. [See the customization guide](customizing-analysis.md).
  * Output files
    * `auspice/<build_name>.json`: output file for visualization in Auspice where `<build_name>` is the name of a build defined in the workflow configuration file.

## Input files

  * `data/metadata.tsv`: tab-delimited description of strain (i.e., sample) attributes
  * `data/sequences.fasta`: genomic sequences whose ids must match the `strain` column in `metadata.tsv`. [See the data preparation guide](data-prep.md).
  * `defaults/include.txt`: list of strain names to _include_ during subsampling and filtering (one strain name per line)
  * `defaults/exclude.txt`: list of strain names to _exclude_ during subsampling and filtering (one strain name per line)

## Output files and directories

  * `auspice/<build_name>.json`: output file for visualization in Auspice where `<build_name>` is the name of your build in the workflow configuration file.
  * `results/aligned.fasta`, `results/filtered.fasta`, etc.: raw results files (dependencies) that are shared across all builds.
  * `results/<build_name>/`: raw results files (dependencies) that are specific to a single build.
  * `logs/`: Log files with error messages and other information about the run.
  * `benchmarks/`: Run-times (and memory usage on Linux systems) for each rule in the workflow.

## Workflow configuration files we might want to customize

  * `my_profiles/<your_profile>/builds.yaml`: workflow configuration file where you can define and configure all the builds you'd like to run. [See the customization guide](customizing-analysis.md).
  * `my_profiles/<your_profile>/config.yaml`: [Snakemake profile configuration](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) where you can define the number of cores to use at once, etc. [See the customization guide](customizing-analysis.md).
  * `defaults/parameters.yaml`: default workflow configuration parameters. Override these settings in the workflow configuration file (`builds.yaml`) above.
  * `defaults/auspice_config.json`: default visualization configuration file. Override these settings in `my_profiles/<your_profile>/auspice_config.json`. [See the customization guide for visualizations](customizing-visualization.md).

## Workflow configuration files we don't need to touch

  * `Snakefile`: entry point for Snakemake commands that also validates inputs. No modification needed.
  * `workflow/snakemake_rules/main_workflow.smk`: defines rules for running each step in the analysis. Modify your `builds.yaml` file, rather than hardcode changes into the snakemake file itself.
  * `workflow/envs/nextstrain.yaml`: specifies computing environment needed to run workflow with the `--use-conda` flag. No modification needed.
  * `workflow/schemas/config.schema.yaml`: defines format (e.g., required fields and types) for  `builds.yaml` files. This can be a useful reference, but no modification needed.
  * `scripts/`: helper scripts for common tasks. No modification needed.
