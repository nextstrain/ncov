rule download_metadata:
    params:
        metadata_url="s3://nextstrain-ncov-private/metadata.tsv.zst",
    output:
        metadata="data/metadata.tsv.zst",
    shell:
        r"""
        aws s3 cp {params.metadata_url} {output.metadata}
        """

rule filter_metadata:
    input:
        metadata="data/metadata.tsv.zst",
        include = config["files"]["include"],
    output:
        metadata="data/prefiltered_metadata.tsv",
    params:
        max_sequences=500000,
        group_by="division year month",
    shell:
        r"""
        augur filter \
            --metadata {input.metadata} \
            --subsample-max-sequences {params.max_sequences} \
            --include {input.include} \
            --group-by {params.group_by} \
            --output-metadata {output.metadata}
        """
