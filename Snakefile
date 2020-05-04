from os import environ
from socket import getfqdn
from getpass import getuser
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# For information on how to run Nextstrain 'regions' runs, see rules/nextstrain_exports.smk
# Regions can be defined as a list or a space-delimited string that is converted
# to a list.
if "regions" not in config:
    REGIONS = ["global"]
elif isinstance(config["regions"], list):
    REGIONS = config["regions"]
else:
    REGIONS = config["regions"].split(" ")

# Define patterns we expect for wildcards.
wildcard_constraints:
    region = "[-a-zA-Z]+",
    date = "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

localrules: download

# Define the format of the output path per region.
# This format can be changed in the future to group builds by different criteria (e.g., division, country, etc.).
REGION_PATH = "results/region/{region}/"

# Create a standard ncov build for auspice, by default.
rule all:
    input:
        auspice_json = expand("auspice/ncov_{region}.json", region=REGIONS),
        tip_frequencies_json = expand("auspice/ncov_{region}_tip-frequencies.json", region=REGIONS)

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

# Include rules to handle primary build logic from multiple sequence alignment
# to output of auspice JSONs for a default build.
include: "rules/builds.smk"

# Include rules specific to the Nextstrain team including custom exports used in
# narratives, etc.
include: "rules/nextstrain_exports.smk"
