"""Small, shared functions used to generate inputs and parameters.
"""
def _get_metadata_by_region(region):
    """Returns a path associated with the metadata for the given region.

    The path can include wildcards that must be provided by the caller through
    the Snakemake `expand` function or through string formatting with `.format`.
    """
    if region == "global":
        return rules.download.output.metadata
    else:
        return rules.adjust_metadata_regions.output.metadata

def _get_metadata_by_wildcards(wildcards):
    """Returns a metadata path based on the given wildcards object.

    This function is designed to be used as an input function.
    """
    return _get_metadata_by_region(wildcards.region)
