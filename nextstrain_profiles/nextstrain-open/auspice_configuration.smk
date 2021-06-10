from snakemake.logging import logger

# Auspice configuration for the (genbank / open) nextstrain build is predominantly through
# a per-build `auspice_config.json` file. The rule defined here creates these files,
# to avoid a plethora of nearly-identical per-build config files stored in the profile.
# NOTE there are two config parameters specified via command-line args to `augur export`:
# (1) The (auspice) title. The main workflow can source this from a variety of locations, but
#     not from within the `auspice_config.json`. For the open runs we expect this to be
#     present in the `config.builds` definitions, but as these are modified by `build_sizes`
#     we create this value here.
# (2) The description (auspice footer). Here the config file specifies one file
#     to be used for all builds.

## Make the title by adding to each build definition, and add an `auspice_config` key
## which specifies a file to be created by the rule defined in this file
if not "build_sizes" in config:
    logger.error("The \"open\" build must have `build_sizes` specified in the `builds.yaml`")
    sys.exit(1)
logger.info("Adding auspice configuration parameters to each build definition")
for build_name, build_info in config["builds"].items():
    pretty_build_name = build_name.replace(f"-{build_info['build_size']['name']}", "").replace("-", " ").title()
    size = build_info['build_size']['n']
    if pretty_build_name == "Global":
        build_info["title"] = f"Genomic epidemiology of novel coronavirus targetting {size} Global sequences"
    else:
        build_info["title"] = f"Genomic epidemiology of novel coronavirus focusing on {pretty_build_name} & targetting {size} sequences"
    build_info['auspice_config'] = f"results/{build_name}/auspice_config.json"

base_config = {
    "title": "This title is not used (but is needed for validation)",
    "build_url": "https://github.com/nextstrain/ncov",
    "maintainers": [
        {"name": "the Nextstrain team", "url": "https://nextstrain.org/"}
    ],
    "data_provenance": [
        {"name": "GENBANK", "url": "https://www.ncbi.nlm.nih.gov/genbank/"}
    ],
    "colorings": [
        { "key": "emerging_lineage", "title": "Emerging Lineage", "type": "categorical" },
        { "key": "clade", "title": "Clade", "type": "categorical"},
        { "key": "S1_mutations", "title": "S1 mutations", "type": "continuous"},
        { "key": "logistic_growth", "title": "Logistic growth", "type": "continuous"},
        { "key": "current_frequency", "title": "Current frequency", "type": "continuous"},
        { "key": "division", "title": "Admin Division", "type": "categorical" },
        { "key": "country", "title": "Country", "type": "categorical" },
        { "key": "region", "title": "Region", "type": "categorical" },
        { "key": "host", "title": "Host", "type": "categorical" },
        { "key": "author", "title": "Authors", "type": "categorical" },
        { "key": "recency", "title": "Submission Date", "type": "categorical" },
        { "key": "genbank_accession", "type": "categorical" }
    ],
    "geo_resolutions": [
        "division",
        "country",
        "region"
    ],
    "display_defaults": {
        "color_by": "clade",
        "distance_measure": "num_date",
        "geo_resolution": "SET_IN_RULE",
        "map_triplicate": "SET_IN_RULE",
        "branch_label": "clade",
        "transmission_lines": False
    },
    "filters": [
        "recency",
        "emerging_lineage",
        "clade",
        "region",
        "country",
        "division",
        "host"
    ],
    "panels": [
        "tree",
        "map",
        "entropy",
        "frequencies"
    ]
}

rule make_auspice_config:
    message: "Making auspice config file {output.auspice_config}"
    output:
        auspice_config = "results/{build_name}/auspice_config.json"
    run:
        import json
        from copy import deepcopy
        config = deepcopy(base_config)
        config["display_defaults"]["geo_resolution"] = "division" if (wildcards.build_name.startswith("north-america") or wildcards.build_name.startswith("oceania")) else "country"
        config["display_defaults"]["map_triplicate"] = wildcards.build_name.startswith("global")
        with open(output.auspice_config, "w") as fh:
            json.dump(config, fh, indent=2)