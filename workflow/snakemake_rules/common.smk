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

def _get_subsampling_scheme_by_build_name(build_name):
    return config["builds"][build_name].get("subsampling_scheme", build_name)

def _get_filter_value(wildcards, key):
    default = config["filter"].get(key, "")
    if wildcards["origin"] == "":
        return default
    return config["filter"].get(wildcards["origin"], {}).get(key, default)

def _get_path_for_input(stage, origin_wildcard):
    """
    A function called to define an input for a Snakemake rule
    This function always returns a local filepath, the format of which lets snakemake decide
    whether to create it (via another rule) or use is as-is.
    """
    path_or_url = config.get("inputs", {}).get(origin_wildcard, {}).get(stage, "")
    scheme = urlsplit(path_or_url).scheme
    remote = bool(scheme)

    # Following checking should be the remit of the rule which downloads the remote resource
    if scheme and scheme!="s3":
        raise Exception(f"Input defined scheme {scheme} which is not yet supported.")

    if stage=="metadata":
        if not path_or_url:
            raise Exception(f"ERROR: config->input->{origin_wildcard}->metadata is not defined.")
        return f"data/downloaded_{origin_wildcard}.tsv" if remote else path_or_url
    if stage=="sequences":
        if not path_or_url:
            raise Exception(f"ERROR: config->input->{origin_wildcard}->sequences is not defined.")
        return f"data/downloaded_{origin_wildcard}.fasta.xz" if remote else path_or_url
    if stage=="aligned":
        if remote:
            return f"results/precomputed-aligned_{origin_wildcard}.fasta"
        elif path_or_url:
            return path_or_url
        else:
            return f"results/aligned_{origin_wildcard}.fasta.xz"
    if stage=="masked":
        if remote:
            return f"results/precomputed-masked_{origin_wildcard}.fasta"
        elif path_or_url:
            return path_or_url
        else:
            return f"results/masked_{origin_wildcard}.fasta.xz"
    if stage=="filtered":
        if remote:
            return f"results/precomputed-filtered_{origin_wildcard}.fasta"
        elif path_or_url:
            return path_or_url
        else:
            return f"results/filtered_{origin_wildcard}.fasta.xz"

    raise Exception(f"_get_path_for_input with unknown stage \"{stage}\"")


def _get_unified_metadata(wildcards):
    """
    Returns a single metadata file representing the input metadata file(s).
    If there was only one supplied metadata file in the `config["inputs"] dict`,
    then that file is returned. Else "results/combined_metadata.tsv" is returned
    which will run the `combine_input_metadata` rule to make it.
    """
    if len(list(config["inputs"].keys()))==1:
        return "results/sanitized_metadata_{origin}.tsv.xz".format(origin=list(config["inputs"].keys())[0])
    return "results/combined_metadata.tsv.xz"

def _get_unified_alignment(wildcards):
    if len(list(config["inputs"].keys()))==1:
        return _get_path_for_input("filtered", list(config["inputs"].keys())[0])
    return "results/combined_sequences_for_subsampling.fasta.xz",

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

def _get_min_date_for_frequencies(wildcards):
    if "frequencies" in config and "min_date" in config["frequencies"]:
        return config["frequencies"]["min_date"]
    else:
        # If not explicitly specified, default to 1 year back from the present
        min_date_cutoff = datetime.date.today() - datetime.timedelta(weeks=52)
        return numeric_date(
            min_date_cutoff
        )

def _get_max_date_for_frequencies(wildcards):
    if "frequencies" in config and "max_date" in config["frequencies"]:
        return config["frequencies"]["max_date"]
    else:
        # Allow users to censor the N most recent days to minimize effects of
        # uneven recent sampling.
        recent_days_to_censor = config.get("frequencies", {}).get("recent_days_to_censor", 0)
        offset = datetime.timedelta(days=recent_days_to_censor)

        return numeric_date(
            date.today() - offset
        )

def _get_upload_inputs(wildcards):
    # The main workflow supports multiple inputs/origins, but our desired file
    # structure under data.nextstrain.org/files/ncov/open/… is designed around
    # a single input/origin.  Intermediates (aligned, masked, filtered, etc)
    # are specific to each input/origin and thus do not match our desired
    # structure, while builds (global, europe, africa, etc) span all
    # inputs/origins (and thus do).  In our desired outcome, the two kinds of
    # files are comingled.  Without changing the main workflow, the mismatch
    # has to be reconciled by the upload rule.  Thus, the upload rule enforces
    # and uses only single-origin configurations.  How third-wave.
    #   -trs, 7 May 2021
    if len(config["S3_DST_ORIGINS"]) != 1:
        raise Exception(f'The "upload" rule requires a single value in S3_DST_ORIGINS (got {config["S3_DST_ORIGINS"]!r}).')

    origin = config["S3_DST_ORIGINS"][0]

    # mapping of remote → local filenames
    uploads = {
        f"aligned.fasta.xz":              f"results/aligned_{origin}.fasta.xz",              # from `rule align`
        f"masked.fasta.xz":               f"results/masked_{origin}.fasta.xz",               # from `rule mask`
        f"filtered.fasta.xz":             f"results/filtered_{origin}.fasta.xz",             # from `rule filter`
        f"mutation-summary.tsv.xz":       f"results/mutation_summary_{origin}.tsv.xz",       # from `rule mutation_summary`
    }

    for build_name in config["builds"]:
        uploads.update({
            f"{build_name}/sequences.fasta.xz": f"results/{build_name}/{build_name}_subsampled_sequences.fasta.xz",   # from `rule combine_samples`
            f"{build_name}/metadata.tsv.xz":    f"results/{build_name}/{build_name}_subsampled_metadata.tsv.xz",      # from `rule combine_samples`
            f"{build_name}/aligned.fasta.xz":   f"results/{build_name}/aligned.fasta.xz",                             # from `rule build_align`
            # export the auspice dataset which matches the subsampled sequences / metadata (see `rule finalize`)
            f"{build_name}/{build_name}.json":                  f"auspice/{config['auspice_json_prefix']}_{build_name}.json",
            f"{build_name}/{build_name}_tip-frequencies.json":  f"auspice/{config['auspice_json_prefix']}_{build_name}_tip-frequencies.json",
            f"{build_name}/{build_name}_root-sequence.json":    f"auspice/{config['auspice_json_prefix']}_{build_name}_root-sequence.json"
        })

    return uploads
