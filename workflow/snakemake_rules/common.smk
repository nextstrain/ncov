"""Small, shared functions used to generate inputs and parameters.
"""
import datetime
from urllib.parse import urlsplit

def numeric_date(dt=None):
    """
    Convert datetime object to the numeric date.
    The numeric date format is YYYY.F, where F is the fraction of the year passed
    Parameters
    ----------
     dt:  datetime.datetime, None
        date of to be converted. if None, assume today
    """
    from calendar import isleap

    if dt is None:
        dt = datetime.datetime.now()

    days_in_year = 366 if isleap(dt.year) else 365
    try:
        res = dt.year + (dt.timetuple().tm_yday-0.5) / days_in_year
    except:
        res = None

    return res

def _trim_origin(origin):
    """the origin wildcard includes a leading `_`. This function returns the value without this `_`"""
    if origin=="":
        return ""
    return origin[1:]

def _get_subsampling_scheme_by_build_name(build_name):
    return config["builds"][build_name].get("subsampling_scheme", build_name)

def _get_filter_value(wildcards, key):
    default = config["filter"].get(key, "")
    if wildcards["origin"] == "":
        return default
    return config["filter"].get(_trim_origin(wildcards["origin"]), {}).get(key, default)

def _get_path_for_input(stage, origin_wildcard):
    """
    A function called to define an input for a Snakemake rule
    This function always returns a local filepath, the format of which decides whether rules should
    create this by downloading from a remote resource, or create it by a local compute rule.
    """
    if not origin_wildcard:
        # No origin wildcards => deprecated single inputs (e.g. `config["sequences"]`) which cannot
        # be downloaded from remote resources
        if config.get("inputs"):
            raise Exception("ERROR: empty origin wildcard but config defines 'inputs`")
        path_or_url = config[stage] if stage in ["metadata", "sequences"] else ""
        remote = False
    else:
        trimmed_origin = _trim_origin(origin_wildcard)
        path_or_url = config.get("inputs", {}).get(trimmed_origin, {}).get(stage, "")
        scheme = urlsplit(path_or_url).scheme
        remote = bool(scheme)

        # Following checking should be the remit of the rule which downloads the remote resource
        if scheme and scheme!="s3":
            raise Exception(f"Input defined scheme {scheme} which is not yet supported.")

        ## Basic checking which could be taken care of by the config schema
        ## If asking for metadata/sequences, the config _must_ supply a `path_or_url`
        if path_or_url=="" and stage in ["metadata", "sequences"]:
            raise Exception(f"ERROR: config->input->{trimmed_origin}->{stage} is not defined.")

    if stage=="metadata":
        return f"data/downloaded{origin_wildcard}.tsv" if remote else path_or_url
    if stage=="sequences":
        return f"data/downloaded{origin_wildcard}.fasta" if remote else path_or_url
    if stage=="prefiltered":
        # we don't expose the option to download a prefiltered alignment - it must be computed locally
        return f"results/prefiltered{origin_wildcard}.fasta"
    if stage=="aligned":
        return f"results/precomputed-aligned{origin_wildcard}.fasta" if remote else f"results/aligned{origin_wildcard}.fasta"
    if stage=="to-exclude":
        return f"results/precomputed-to-exclude{origin_wildcard}.txt" if remote else f"results/to-exclude{origin_wildcard}.txt"
    if stage=="aligned-filtered": # refiltering. Different to filtered!
        return f"results/precomputed-aligned-filtered{origin_wildcard}.fasta" if remote else f"results/aligned-filtered{origin_wildcard}.fasta"
    if stage=="masked":
        return f"results/precomputed-masked{origin_wildcard}.fasta" if remote else f"results/masked{origin_wildcard}.fasta"
    if stage=="filtered":
        return f"results/precomputed-filtered{origin_wildcard}.fasta" if remote else f"results/filtered{origin_wildcard}.fasta"

    raise Exception(f"_get_path_for_input with unknown stage \"{stage}\"")


def _get_unified_metadata(wildcards):
    """
    Returns a single metadata file representing the input metadata file(s).
    If there was only one supplied metadata file (e.g. the deprecated
    `config["metadata"]` syntax, or one entry in the `config["inputs"] dict`)
    then that file is returned. Else "results/combined_metadata.tsv" is returned
    which will run the `combine_input_metadata` rule to make it.
    """
    if not config.get("inputs"):
        return config["metadata"]
    if len(list(config["inputs"].keys()))==1:
        return _get_path_for_input("metadata", "_"+list(config["inputs"].keys())[0])
    return "results/combined_metadata.tsv"

def _get_unified_alignment(wildcards):
    if not config.get("inputs"):
        return "results/filtered.fasta"
    if len(list(config["inputs"].keys()))==1:
        return _get_path_for_input("filtered", "_"+list(config["inputs"].keys())[0])
    return "results/combined_sequences_for_subsampling.fasta",

def _get_metadata_by_build_name(build_name):
    """Returns a path associated with the metadata for the given build name.

    The path can include wildcards that must be provided by the caller through
    the Snakemake `expand` function or through string formatting with `.format`.
    """
    if build_name == "global" or "region" not in config["builds"][build_name]:
        return _get_unified_metadata({})
    else:
        return rules.adjust_metadata_regions.output.metadata

def _get_metadata_by_wildcards(wildcards):
    """Returns a metadata path based on the given wildcards object.

    This function is designed to be used as an input function.
    """
    return _get_metadata_by_build_name(wildcards.build_name)

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

def _get_max_date_for_frequencies(wildcards):
    if "frequencies" in config and "max_date" in config["frequencies"]:
        return config["frequencies"]["max_date"]
    else:
        return numeric_date(date.today())

def _get_first(config, *keys):
    """
    Get the value of the first key in *keys* that exists in *config* and has a
    non-empty value.

    Raises a :class:`KeyError` if none of the *keys* exist with a suitable
    value.
    """
    for key in keys:
        if config.get(key) not in {"", None}:
            return config[key]
    raise KeyError(str(keys))
