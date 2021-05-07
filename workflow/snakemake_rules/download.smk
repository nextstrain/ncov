# Here we define a number of rules which produce the same output as rules in
# the main_workflow, but do so by downloading the files from S3 rather than
# computing them.

# These rules are run by input functions requesting their filenames, which differ
# from those which would be produced locally.
# For instance, the following rule:
#       rule: align
#           input: _get_path_for_input
#           ...
# will result in an input file looking like "results/aligned_{origin}.fasta" or
# "results/download-aligned_{origin}.fasta" (which one is chosen depends on the
# supplied `config`). In the latter case, `rule download_aligned` will be used.
# See https://github.com/nextstrain/ncov/compare/remote-files for an example of
# how we could leverage snakemake to do this without needing a separate rule!


def _infer_decompression(input):
    """
    Returns a shell command to decompress the piped stream,
    which will itself produce a stream of decompressed data to stdout.
    If no decompression is needed, returns `cat`.
    NOTE: a lot of this will become unnecessary once `augur` handles
    compressed sequence inputs.
    """
    if input.endswith(".xz"):
        return "xz -dcq"
    if input.endswith(".gz"):
        return "gunzip -cq"
    return "cat"


rule download_sequences:
    message: "Downloading sequences from {params.address} -> {output.sequences}"
    output:
        sequences = "data/downloaded_{origin}.fasta"
    conda: config["conda_environment"]
    params:
        address = lambda w: config["inputs"][w.origin]["sequences"],
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["sequences"])
    shell:
        """
        aws s3 cp {params.address} - | {params.deflate} > {output.sequences:q}
        """

rule download_metadata:
    message: "Downloading metadata from {params.address} -> {output.metadata}"
    output:
        metadata = "data/downloaded_{origin}.tsv"
    conda: config["conda_environment"]
    params:
        address = lambda w: config["inputs"][w.origin]["metadata"],
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["metadata"]),
    shell:
        """
        aws s3 cp {params.address} - | {params.deflate} > {output.metadata:q}
        """

rule download_aligned:
    message: "Downloading aligned fasta files from {params.address} -> {output.sequences}"
    output:
        sequences = "results/precomputed-aligned_{origin}.fasta"
    conda: config["conda_environment"]
    params:
        address = lambda w: config["inputs"][w.origin]["aligned"],
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["aligned"])
    shell:
        """
        aws s3 cp {params.address} - | {params.deflate} > {output.sequences:q}
        """

rule download_diagnostic:
    message: """
    Downloading diagnostic files:
        {params.diagnostics_address} -> {output.diagnostics}
        {params.flagged_address} -> {output.flagged}
        {params.to_exclude_address} -> {output.to_exclude}
    """
    output:
        diagnostics = "results/precomputed-sequence-diagnostics_{origin}.tsv",
        flagged = "results/precomputed-flagged-sequences_{origin}.tsv",
        to_exclude = "results/precomputed-to-exclude_{origin}.txt"
    conda: config["conda_environment"]
    params:
        # Only `to-exclude` is defined via the config, so we make some assumptions about the format of the other filenames
        to_exclude_address = lambda w: config["inputs"][w.origin]["to-exclude"],
        flagged_address = lambda w: config["inputs"][w.origin]["to-exclude"].replace(f'to-exclude{w.origin}.txt', f'flagged-sequences{w.origin}.tsv'),
        diagnostics_address = lambda w: config["inputs"][w.origin]["to-exclude"].replace(f'to-exclude{w.origin}.txt', f'sequence-diagnostics{w.origin}.tsv'),
        # assume the compression is the same across all 3 addresses
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["to-exclude"])
    shell:
        """
        aws s3 cp {params.to_exclude_address} - | {params.deflate} > {output.to_exclude:q}
        aws s3 cp {params.flagged_address} - | {params.deflate} > {output.flagged:q}
        aws s3 cp {params.diagnostics_address} - | {params.deflate} > {output.diagnostics:q}
        """


rule download_masked:
    message: "Downloading aligned & masked FASTA from {params.address} -> {output.sequences}"
    output:
        sequences = "results/precomputed-masked_{origin}.fasta"
    conda: config["conda_environment"]
    params:
        address = lambda w: config["inputs"][w.origin]["masked"],
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["masked"])
    shell:
        """
        aws s3 cp {params.address} - | {params.deflate} > {output.sequences:q}
        """


rule download_filtered:
    message: "Downloading pre-computed filtered alignment from {params.address} -> {output.sequences}"
    output:
        sequences = "results/precomputed-filtered_{origin}.fasta"
    conda: config["conda_environment"]
    params:
        address = lambda w: config["inputs"][w.origin]["filtered"],
        deflate = lambda w: _infer_decompression(config["inputs"][w.origin]["filtered"])
    shell:
        """
        aws s3 cp {params.address} - | {params.deflate} > {output.sequences:q}
        """
