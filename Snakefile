from getpass import getuser
from shutil import which
from socket import getfqdn
from snakemake.logging import logger
from snakemake.utils import validate
import sys

#
# Verify that required versions of dependencies are installed.
#
if not (which("docker") and which("nextstrain")) and not which("conda"):
    logger.error("ERROR: The Nextstrain CLI and Docker must be installed or conda must be installed.")
    logger.error("Please follow installation instructions at https://nextstrain.org/docs/ and try again.")
    sys.exit(1)

# Load and validate configuration.
configfile: "config/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# For information on how to run 'regions' runs, see rules/regions.smk
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

include: "rules/utils.smk"
include: "rules/builds.smk"
include: "rules/regions.smk"
