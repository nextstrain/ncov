rule download_metadata:
    params:
        metadata_url="s3://nextstrain-ncov-private/metadata.tsv.zst",
    output:
        metadata="data/metadata.tsv.zst",
    shell:
        """
        aws s3 cp {params.metadata_url} {output.metadata}
        """

rule filter_metadata:
    input:
        metadata="data/metadata.tsv.zst",
    output:
        metadata="data/prefiltered_metadata.tsv",
    params:
        max_sequences=500000,
        group_by="division year month",
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --subsample-max-sequences {params.max_sequences} \
            --group-by {params.group_by} \
            --output-metadata {output.metadata}
        """
