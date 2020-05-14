from os import environ
from socket import getfqdn
from getpass import getuser
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

LOCATIONS = [
    (location_type, location)
    for location_type, location_values in config["locations"].items()
    for location in location_values
]
LOCATION_TYPES, LOCATION_NAMES = list(zip(*LOCATIONS))

# Define patterns we expect for wildcards.
wildcard_constraints:
    location_type = "[-a-zA-Z]+",
    location_name = "[-a-zA-Z]+",
    date = "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"

localrules: download

# Create a standard ncov build for auspice, by default.
rule all:
    input:
        auspice_json = expand("auspice/ncov_{location_type}_{location_name}.json", zip, location_type=LOCATION_TYPES, location_name=LOCATION_NAMES),
        tip_frequencies_json = expand("auspice/ncov_{location_type}_{location_name}_tip-frequencies.json", zip, location_type=LOCATION_TYPES, location_name=LOCATION_NAMES)

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
