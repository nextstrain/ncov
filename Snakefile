from os import environ
from socket import getfqdn
from getpass import getuser
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# In our shell rules, we try to use the same Python executable running
# Snakemake, if known, by interpolating {python:q}.  The aim is to prevent
# issues with multiple Python versions or symlinks for python vs. python3.
python = sys.executable or "python3"

# default build if none specified in config
if "builds" not in config:
    config["builds"] = {
        "global":{
                    "geographic_scale": "region",
                    "subsampling_scheme": "region_global",
                 }
        }

BUILD_NAMES = list(config["builds"].keys())

# Define patterns we expect for wildcards.
wildcard_constraints:
    # Allow build names to contain alpha characters, underscores, and hyphens
    # but not special strings used for Nextstrain builds.
    build_name = r'(?:[_a-zA-Z-](?!(tip-frequencies|gisaid|zh)))+',
    date = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

localrules: download

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

# Include small, shared functions that help build inputs and parameters.
include: "rules/common.smk"

# Include rules to handle primary build logic from multiple sequence alignment
# to output of auspice JSONs for a default build.
include: "rules/builds.smk"

# Include rules specific to the Nextstrain team including custom exports used in
# narratives, etc.
include: "rules/nextstrain_exports.smk"

# Include a custom Snakefile that specifies `localrules` required by the user's
# workflow environment.
if "localrules" in config:
    include: config["localrules"]
