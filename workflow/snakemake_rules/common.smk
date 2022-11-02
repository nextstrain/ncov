"""Small, shared functions used to generate inputs and parameters.
"""
import datetime
from itertools import product
from shlex import (
    quote as shquote,       # shquote() is used in this file and also other workflow files
    split as shsplitwords,
)
from urllib.parse import urlsplit

def shquotewords(s: str) -> str:
    """
    Split string *s* into (POSIX) shell words, quote each word, and join them
    back into a string.

    This is suitable for properly quoting multi-word, user-defined values which
    should follow shell quoting and escaping semantics (e.g. to allow spaces in
    single words) but not allow shell features like variable interpolation,
    command substition, redirection, piping, etc.
    """
    return " ".join(shquote(word) for word in shsplitwords(s))

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
    return config["builds"].get(build_name, {}).get("subsampling_scheme", build_name)

def _get_skipped_inputs_for_diagnostic(wildcards):
    """Build an argument for the diagnostic script with a list of inputs to skip.
    """
    inputs = config["inputs"]
    diagnostics_key = "skip_diagnostics"

    arg_parts = []
    for input_name in inputs.keys():
        skip_diagnostics = config["filter"].get(diagnostics_key, False)

        if input_name in config["filter"] and diagnostics_key in config["filter"][input_name]:
            skip_diagnostics = config["filter"][input_name][diagnostics_key]

        if skip_diagnostics:
            arg_parts.append(input_name)

    if len(arg_parts) > 0:
        argument = f"--skip-inputs {' '.join(arg_parts)}"
    else:
        argument = ""

    return argument

def _get_filter_min_length_query(wildcards):
    """Build a sequence length filter query for each input, checking for
    input-specific length requirements.
    """
    inputs = config["inputs"]
    length_key = "min_length"

    query_parts = []
    for input_name in inputs.keys():
        min_length = config["filter"][length_key]

        if input_name in config["filter"] and length_key in config["filter"][input_name]:
            min_length = config["filter"][input_name][length_key]

        # We only annotate input-specific columns on metadata when there are
        # multiple inputs.
        if len(inputs) > 1:
            # Input names can contain characters that make them invalid Python
            # variable names. As such, we escape the column names with backticks
            # as recommended by pandas:
            #
            # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html
            query_parts.append(f"(`{input_name}` == 'yes' & _length >= {min_length})")
        else:
            query_parts.append(f"(_length >= {min_length})")

    query = " | ".join(query_parts)
    return f"--query {shquote(query)}"

def _get_filter_value(wildcards, key):
    for input_name in config["inputs"].keys():
        if input_name in config["filter"] and key in config["filter"][input_name]:
            print(
                f"ERROR: We no longer support input-specific filtering with the '{key}' parameter.",
                "Remove this parameter from your configuration file and try running the workflow again.",
                file=sys.stderr,
            )
            sys.exit(1)

    return config["filter"].get(key, "")

def _get_path_for_input(stage, origin_wildcard):
    """
    A function called to define an input for a Snakemake rule
    This function always returns a local filepath, the format of which lets snakemake decide
    whether to create it (via another rule) or use is as-is.
    """
    input_file = config.get("inputs", {}).get(origin_wildcard, {}).get(stage, "")

    if input_file:
        return path_or_url(input_file, keep_local=True)

    if stage in {"metadata", "sequences"}:
        raise Exception(f"ERROR: config->input->{origin_wildcard}->{stage} is not defined.")
    elif stage in {"aligned"}:
        return f"results/{stage}_{origin_wildcard}.fasta.xz"
    else:
        raise Exception(f"_get_path_for_input with unknown stage \"{stage}\"")


def _get_unified_metadata(wildcards):
    """
    Returns a single metadata file representing the input metadata file(s).
    If there was only one supplied metadata file in the `config["inputs"] dict`,
    then that file is run through `sanitize_metadata` and the new file name returned.
    Else "results/combined_metadata.tsv.xz" is returned which will run the
    `combine_input_metadata` rule (and `sanitize_metadata` rule) to make it.
    """
    if len(list(config["inputs"].keys()))==1:
        input_name, input_record = list(config["inputs"].items())[0]
        if input_record.get("skip_sanitize_metadata"):
            return _get_path_for_input("metadata", input_name)
        else:
            return "results/sanitized_metadata_{origin}.tsv.xz".format(origin=input_name)

    return "results/combined_metadata.tsv.xz"

def _get_unified_alignment(wildcards):
    if len(list(config["inputs"].keys()))==1:
        return _get_path_for_input("aligned", list(config["inputs"].keys())[0])
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
    if wildcards.build_name in config["frequencies"] and "min_date" in config["frequencies"][wildcards.build_name]:
        return config["frequencies"][wildcards.build_name]["min_date"]
    elif "frequencies" in config and "min_date" in config["frequencies"]:
        return config["frequencies"]["min_date"]
    else:
        # If not explicitly specified, default to 1 year back from the present
        min_date_cutoff = datetime.date.today() - datetime.timedelta(weeks=52)
        return numeric_date(
            min_date_cutoff
        )

def _get_max_date_for_frequencies(wildcards):
    if wildcards.build_name in config["frequencies"] and "max_date" in config["frequencies"][wildcards.build_name]:
        return config["frequencies"][wildcards.build_name]["max_date"]
    elif "frequencies" in config and "max_date" in config["frequencies"]:
        return config["frequencies"]["max_date"]
    else:
        # Allow users to censor the N most recent days to minimize effects of
        # uneven recent sampling.
        recent_days_to_censor = config.get("frequencies", {}).get("recent_days_to_censor", 0)
        offset = datetime.timedelta(days=recent_days_to_censor)

        return numeric_date(
            datetime.date.today() - offset
        )

def _get_upload_inputs(wildcards):
    # The main workflow supports multiple inputs/origins, but our desired file
    # structure under data.nextstrain.org/files/ncov/open/… is designed around
    # a single input/origin.  Intermediates (aligned, etc)
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

    # This function bakes in these assumptions here about the build names used
    # for the nextstrain.org/ncov/gisaid and …/open builds and then
    # special-cases them below.
    regions = {"global", "africa", "asia", "europe", "north-america", "oceania", "south-america"}
    timespans = {"1m", "2m", "6m", "all-time"}
    region_timespan_builds = [f"{region}_{timespan}" for region, timespan in product(regions, timespans)]

    # mapping of remote → local filenames
    build_files = {}
    for build_name in config["builds"]:
        if build_name in region_timespan_builds:
            region, timespan = build_name.split("_")

            # We name remote files only by region (for now), so only include
            # the 6m timespan builds.
            if timespan != "6m":
                continue

            upload_name = region
        else:
            upload_name = build_name

        build_files.update({
            f"{upload_name}/sequences.fasta.xz": f"results/{build_name}/{build_name}_subsampled_sequences.fasta.xz",   # from `rule combine_samples`
            f"{upload_name}/metadata.tsv.xz":    f"results/{build_name}/{build_name}_subsampled_metadata.tsv.xz",      # from `rule combine_samples`
            f"{upload_name}/aligned.fasta.xz":   f"results/{build_name}/aligned.fasta.xz",                             # from `rule build_align`
            # export the auspice dataset which matches the subsampled sequences / metadata (see `rule finalize`)
            f"{upload_name}/{upload_name}.json":                  f"auspice/{config['auspice_json_prefix']}_{build_name}.json",
            f"{upload_name}/{upload_name}_tip-frequencies.json":  f"auspice/{config['auspice_json_prefix']}_{build_name}_tip-frequencies.json",
            f"{upload_name}/{upload_name}_root-sequence.json":    f"auspice/{config['auspice_json_prefix']}_{build_name}_root-sequence.json"
        })
    return build_files
