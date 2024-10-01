rule clades_21L:
    input:
        clades = "defaults/clades.tsv",
        exclude_clades = "nextstrain_profiles/nextstrain-gisaid-21L/exclude-clades.tsv",
    output:
        clades = "results/clades_21L.tsv",
    log: "logs/clades_21L.txt"
    benchmark: "benchmarks/clades_21L.txt"
    conda: config["conda_environment"]
    shell:
        r"""
        exec 2> {log:q}

          ./scripts/expand-clade-definitions {input.clades:q} \
        | tsv-join \
            --header \
            --exclude \
            --filter-file {input.exclude_clades:q} \
            --key-fields clade \
        > {output.clades:q}
        """


rule gisaid_21L_metadata:
    input:
        references = "data/references_metadata.tsv",
        metadata = path_or_url("s3://nextstrain-ncov-private/metadata.tsv.zst", keep_local=True),
        exclude_clades = "nextstrain_profiles/nextstrain-gisaid-21L/exclude-clades.tsv",
    output:
        metadata = "results/gisaid_21L_metadata.tsv.zst",
    log: "logs/gisaid_21L_metadata.txt"
    benchmark: "benchmarks/gisaid_21L_metadata.txt"
    conda: config["conda_environment"]
    threads: 8
    shell:
        r"""
        exec 2> {log:q}

        ./scripts/tsv-cast-header \
            <(unzstd < {input.metadata:q}) \
            {input.references:q} \
        | zstd \
        > {output.metadata:q}

        < {input.metadata:q} \
          unzstd \
        | tsv-join \
            --header \
            --exclude \
            --filter-file {input.exclude_clades:q} \
            --key-fields clade \
            --data-fields Nextstrain_clade \
        | sed 1d \
        | zstd -T$(({threads} - 2)) \
        >> {output.metadata:q}
        """


rule gisaid_21L_strains:
    input:
        metadata = "results/gisaid_21L_metadata.tsv.zst",
    output:
        strains = "results/gisaid_21L_strains.txt",
    log: "logs/gisaid_21L_strains.txt"
    benchmark: "benchmarks/gisaid_21L_strains.txt"
    conda: config["conda_environment"]
    shell:
        r"""
        exec 2> {log:q}

        < {input.metadata:q} \
          unzstd \
        | tsv-select --header -f strain \
        | sed 1d \
        > {output.strains:q}
        """


rule gisaid_21L_aligned:
    input:
        references = "data/references_sequences.fasta",
        aligned = path_or_url("s3://nextstrain-ncov-private/aligned.fasta.zst", keep_local=True),
        strains = "results/gisaid_21L_strains.txt",
    output:
        aligned = "results/gisaid_21L_aligned.fasta.zst",
    log: "logs/gisaid_21L_aligned.txt"
    benchmark: "benchmarks/gisaid_21L_aligned.txt"
    conda: config["conda_environment"]
    threads: 8
    shell:
        r"""
        exec 2> {log:q}

        < {input.references:q} \
          seqkit grep --by-name --pattern 21L \
        | zstd \
        > {output.aligned}

        < {input.aligned:q} \
          unzstd \
        | seqkit grep --by-name -f {input.strains:q} \
        | zstd -T$(({threads} - 2)) \
        >> {output.aligned:q}
        """
