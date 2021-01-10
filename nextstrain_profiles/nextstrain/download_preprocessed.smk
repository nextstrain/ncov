# Here we define a number of rules which produce the same output as rules in
# the main_workflow, but do so by downloading the files from S3 rather than
# computing them.

# We control which rule runs via `ruleorder` declarations
ruleorder: download_aligned > align
ruleorder: download_filtered > filter
ruleorder: download_refiltered > refilter
ruleorder: download_masked > mask
ruleorder: download_diagnostic > diagnostic

rule download_aligned:
    message: "Downloading aligned fasta files from S3 bucket {params.s3_bucket}"
    output:
        sequences = "results/aligned.fasta"
    conda: config["conda_environment"]
    params:
        compression = config['preprocess']['compression'],
        deflate = config['preprocess']['deflate'],
        s3_bucket = config["S3_BUCKET"]
    shell:
        """
        aws s3 cp s3://{params.s3_bucket}/aligned.fasta.{params.compression} - | {params.deflate} > {output.sequences:q}
        """

rule download_diagnostic:
    message: "Downloading diagnostic files from S3 bucket {params.s3_bucket}"
    output:
        diagnostics = "results/sequence-diagnostics.tsv",
        flagged = "results/flagged-sequences.tsv",
        to_exclude = "results/to-exclude.txt"
    conda: config["conda_environment"]
    params:
        compression = config['preprocess']['compression'],
        deflate = config['preprocess']['deflate'],
        s3_bucket = config["S3_BUCKET"]
    shell:
        """
        aws s3 cp s3://{params.s3_bucket}/sequence-diagnostics.tsv.{params.compression} - | {params.deflate} > {output.diagnostics:q}
        aws s3 cp s3://{params.s3_bucket}/flagged-sequences.tsv.{params.compression} - | {params.deflate} > {output.flagged:q}
        aws s3 cp s3://{params.s3_bucket}/to-exclude.txt.{params.compression} - | {params.deflate} > {output.to_exclude:q}
        """

rule download_refiltered:
    message: "Downloading quality filtered files from S3 bucket {params.s3_bucket}"
    output:
        sequences = "results/aligned-filtered.fasta"
    conda: config["conda_environment"]
    params:
        compression = config['preprocess']['compression'],
        deflate = config['preprocess']['deflate'],
        s3_bucket = config["S3_BUCKET"]
    shell:
        """
        aws s3 cp s3://{params.s3_bucket}/aligned-filtered.fasta.{params.compression} - | {params.deflate} > {output.sequences:q}
        """


rule download_masked:
    message: "Downloading aligned masked fasta files from S3 bucket {params.s3_bucket}"
    output:
        sequences = "results/masked.fasta"
    conda: config["conda_environment"]
    params:
        compression = config['preprocess']['compression'],
        deflate = config['preprocess']['deflate'],
        s3_bucket = config["S3_BUCKET"]
    shell:
        """
        aws s3 cp s3://{params.s3_bucket}/masked.fasta.{params.compression} - | {params.deflate} > {output.sequences:q}
        """


rule download_filtered:
    message: "Downloading final filtered fasta files from S3 bucket {params.s3_bucket}"
    output:
        sequences = "results/filtered.fasta"
    conda: config["conda_environment"]
    params:
        compression = config['preprocess']['compression'],
        deflate = config['preprocess']['deflate'],
        s3_bucket = config["S3_BUCKET"]
    shell:
        """
        aws s3 cp s3://{params.s3_bucket}/filtered.fasta.{params.compression} - | {params.deflate} > {output.sequences:q}
        """