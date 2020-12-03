"""Small, shared functions used to generate inputs and parameters.
"""
def _get_subsampling_scheme_by_build_name(build_name):
    return config["builds"][build_name].get("subsampling_scheme", build_name)

def _get_main_metadata(wildcards):
    """Returns the appropriate metadata file for all analysis steps"""
    if isinstance(config["metadata"], str):
        return config["metadata"]
    else:
        return "results/combined_metadata.tsv"


def _get_metadata_by_build_name(build_name):
    """Returns a path associated with the metadata for the given build name.

    The path can include wildcards that must be provided by the caller through
    the Snakemake `expand` function or through string formatting with `.format`.
    """
    if build_name == "global" or "region" not in config["builds"][build_name]:
        return _get_main_metadata({})
    else:
        return rules.adjust_metadata_regions.output.metadata

def _get_sequences_by_wildcards(wildcards):
    if isinstance(config["sequences"], str):
        return config["sequences"]
    return config["sequences"][wildcards['origin']]

def _get_metadata_by_wildcards(wildcards):
    """Returns a metadata path based on the given wildcards object.

    This function is designed to be used as an input function.
    """
    # When run for rules involved in subsampled build creation, this returns the subsampled metadata
    if "build_name" in wildcards.keys():
        return _get_metadata_by_build_name(wildcards.build_name)
    # When run for rules involving inital aligning/filtering steps, this returns the correct input metatata file
    elif "origin" in wildcards.keys():
        if isinstance(config["metadata"], str):
            return config["metadata"]
        return config["metadata"][wildcards['origin']]
    raise Exception("_get_metadata_by_wildcards called with unknown wildcards!")

def _get_filter_param(setting):
    """get filter param which may be specific to the origin wildcard"""
    def _get_param(wildcards):
        if setting in config['filter'].get(wildcards.origin, {}):
            return config['filter'][wildcards.origin][setting]
        return config['filter'][setting]
    return _get_param

def _get_sampling_trait_for_wildcards(wildcards):
    if wildcards.build_name in config["exposure"]:
        return config["exposure"][wildcards.build_name]["trait"]
    else:
        return config["exposure"]["default"]["trait"]

def _get_exposure_trait_for_wildcards(wildcards):
    if wildcards.build_name in config["exposure"]:
        return config["exposure"][wildcards.build_name]["exposure"]
    else:
        return config["exposure"]["default"]["exposure"]

def _get_trait_columns_by_wildcards(wildcards):
    if wildcards.build_name in config["traits"]:
        return config["traits"][wildcards.build_name]["columns"]
    else:
        return config["traits"]["default"]["columns"]

def _get_sampling_bias_correction_for_wildcards(wildcards):
    if wildcards.build_name in config["traits"] and 'sampling_bias_correction' in config["traits"][wildcards.build_name]:
        return config["traits"][wildcards.build_name]["sampling_bias_correction"]
    else:
        return config["traits"]["default"]["sampling_bias_correction"]
