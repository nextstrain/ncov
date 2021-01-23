import copy
from datetime import date
import os
from os import environ
from socket import getfqdn
from getpass import getuser
from snakemake.logging import logger
from snakemake.utils import validate
import time

# Store the user's configuration prior to loading defaults, so we can check for
# reused subsampling scheme names in the user's config. We need to make a deep
# copy because Snakemake will deep merge the subsampling dictionary later,
# modifying the values of a reference or shallow copy. Note that this loading of
# the user's config prior to the defaults below depends on the order Snakemake
# loads its configfiles. Specifically, the order of config loading is:
#
# 1. First, configfile arguments are loaded and config is built from these [1].
# 2. Then, config arguments are loaded and override existing config keys [2].
# 3. Then, the Snakefile is parsed and configfile directive inside the Snakefile is processed [3].
#    When configfile is loaded from the directive in the Snakefile, the config
#    dictionary is deep merged with the files [4] from the externally provided
#    config files. This is the only place the deep merge happens using the
#    update_config function [5].
#
# [1] https://github.com/snakemake/snakemake/blob/a7ac40c96d6e2af47102563d0478a2220e2a2ab7/snakemake/__init__.py#L466-L471
# [2] https://github.com/snakemake/snakemake/blob/a7ac40c96d6e2af47102563d0478a2220e2a2ab7/snakemake/__init__.py#L476-L477
# [3] https://github.com/snakemake/snakemake/blob/a7ac40c96d6e2af47102563d0478a2220e2a2ab7/snakemake/__init__.py#L551-L553
# [4] https://github.com/snakemake/snakemake/blob/a7ac40c96d6e2af47102563d0478a2220e2a2ab7/snakemake/workflow.py#L1088-L1094
# [5] https://github.com/snakemake/snakemake/blob/a7ac40c96d6e2af47102563d0478a2220e2a2ab7/snakemake/utils.py#L455-L476
user_subsampling = copy.deepcopy(config.get("subsampling", {}))

configfile: "defaults/parameters.yaml"
validate(config, schema="workflow/schemas/config.schema.yaml")

# Check for overlapping subsampling schemes in user and default
# configurations. For now, issue a deprecation warning, so users know they
# should rename their subsampling schemes. In the future, this reuse of the same
# name will cause an error.
subsampling_config = config.get("subsampling", {})
overlapping_schemes = []
for scheme_name, scheme in user_subsampling.items():
    if scheme_name in subsampling_config and subsampling_config.get(scheme_name) != scheme:
        overlapping_schemes.append(scheme_name)

if len(overlapping_schemes) > 0:
    logger.warning(f"WARNING: The following subsampling scheme(s) have the same name as a default scheme in this workflow but different definitions:")
    logger.warning("")
    for scheme in overlapping_schemes:
        logger.warning(f"  - {scheme}")
    logger.warning("")
    logger.warning("  This means Snakemake will merge your scheme with the default scheme and may produce unexpected behavior.")
    logger.warning(f"  To avoid errors in your workflow, rename your schemes with unique names (e.g., 'custom_{overlapping_schemes[0]}')")
    logger.warning("  In future versions of this workflow, overlapping subsampling scheme names will produce an error.")
    logger.warning("")
    time.sleep(5)

# default build if none specified in config
if "builds" not in config:
    config["builds"] = {
        "global": {
            "subsampling_scheme": "region_global",
        }
    }

# Allow users to specify a list of active builds from the command line.
if config.get("active_builds"):
    BUILD_NAMES = config["active_builds"].split(",")
else:
    BUILD_NAMES = list(config["builds"].keys())

# Construct the correct absolute path to the conda environment based on the
# top-level Snakefile's directory and a path defined in the Snakemake config
# that is relative to that top-level directory.
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
CONDA_ENV_PATH = os.path.join(SNAKEMAKE_DIR, config["conda_environment"])
config["conda_environment"] = CONDA_ENV_PATH

# Define patterns we expect for wildcards.
wildcard_constraints:
    # Allow build names to contain alpha characters, underscores, and hyphens
    # but not special strings used for Nextstrain builds.
    build_name = r'(?:[_a-zA-Z-](?!(tip-frequencies)))+',
    date = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

localrules: download_metadata, download_sequences, download, upload, clean

# Create a standard ncov build for auspice, by default.
rule all:
    input:
        auspice_json = expand("auspice/ncov_{build_name}.json", build_name=BUILD_NAMES),
        tip_frequency_json = expand("auspice/ncov_{build_name}_tip-frequencies.json", build_name=BUILD_NAMES)

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

rule dump_config:
    run:
        import yaml, sys
        yaml.dump(config, sys.stdout, explicit_start = True, explicit_end = True)

# Include small, shared functions that help build inputs and parameters.
include: "workflow/snakemake_rules/common.smk"

# Include rules to handle primary build logic from multiple sequence alignment
# to output of auspice JSONs for a default build.
include: "workflow/snakemake_rules/main_workflow.smk"

# Include a custom Snakefile that specifies `localrules` required by the user's
# workflow environment.
if "localrules" in config:
    include: config["localrules"]

# Include custom rules defined in the config.
if "custom_rules" in config:
    for rule_file in config["custom_rules"]:
        include: rule_file
