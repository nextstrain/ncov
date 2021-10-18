rule sanitize_metadata:
    input:
        metadata=lambda wildcards: _get_path_for_input("metadata", wildcards.origin)
    output:
        metadata="results/sanitized_metadata_{origin}.tsv.xz"
    benchmark:
        "benchmarks/sanitize_metadata_{origin}.txt"
    conda:
        config["conda_environment"]
    log:
        "logs/sanitize_metadata_{origin}.txt"
    params:
        metadata_id_columns=config["sanitize_metadata"]["metadata_id_columns"],
        database_id_columns=config["sanitize_metadata"]["database_id_columns"],
        parse_location_field=f"--parse-location-field {config['sanitize_metadata']['parse_location_field']}" if config["sanitize_metadata"].get("parse_location_field") else "",
        rename_fields=config["sanitize_metadata"]["rename_fields"],
        strain_prefixes=config["strip_strain_prefixes"],
        error_on_duplicate_strains="--error-on-duplicate-strains" if config["sanitize_metadata"].get("error_on_duplicate_strains") else "",
    resources:
        mem_mb=2000
    shell:
        """
        python3 scripts/sanitize_metadata.py \
            --metadata {input.metadata} \
            --metadata-id-columns {params.metadata_id_columns:q} \
            --database-id-columns {params.database_id_columns:q} \
            {params.parse_location_field} \
            --rename-fields {params.rename_fields:q} \
            --strip-prefixes {params.strain_prefixes:q} \
            {params.error_on_duplicate_strains} \
            --output {output.metadata} 2>&1 | tee {log}
        """


rule combine_input_metadata:
    # this rule is intended to be run _only_ if we have defined multiple inputs ("origins")
    message:
        """
        Combining metadata files {input.metadata} -> {output.metadata} and adding columns to represent origin
        """
    input:
        metadata=expand("results/sanitized_metadata_{origin}.tsv.xz", origin=config.get("inputs")),
    output:
        metadata = "results/combined_metadata.tsv.xz"
    params:
        origins = lambda wildcards: list(config["inputs"].keys())
    log:
        "logs/combine_input_metadata.txt"
    benchmark:
        "benchmarks/combine_input_metadata.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine_metadata.py --metadata {input.metadata} --origins {params.origins} --output {output.metadata} 2>&1 | tee {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = lambda wildcards: _get_path_for_input("sequences", wildcards.origin),
        genemap = config["files"]["annotation"],
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "results/aligned_{origin}.fasta.xz",
        insertions = "results/insertions_{origin}.tsv",
        translations = expand("results/translations/seqs_{{origin}}.gene.{gene}.fasta.xz", gene=config.get('genes', ['S']))
    params:
        outdir = "results/translations",
        genes = ','.join(config.get('genes', ['S'])),
        basename = "seqs_{origin}",
        strain_prefixes=config["strip_strain_prefixes"],
        # Strip the compression suffix for the intermediate output from the aligner.
        uncompressed_alignment=lambda wildcards, output: Path(output.alignment).with_suffix(""),
        sanitize_log="logs/sanitize_sequences_{origin}.txt"
    log:
        "logs/align_{origin}.txt"
    benchmark:
        "benchmarks/align_{origin}.txt"
    conda: config["conda_environment"]
    threads: 8
    resources:
        mem_mb=3000
    shell:
        """
        python3 scripts/sanitize_sequences.py \
            --sequences {input.sequences} \
            --strip-prefixes {params.strain_prefixes:q} \
            --output /dev/stdout 2> {params.sanitize_log} \
            | nextalign \
            --jobs={threads} \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences /dev/stdin \
            --output-dir {params.outdir} \
            --output-basename {params.basename} \
            --output-fasta {params.uncompressed_alignment} \
            --output-insertions {output.insertions} > {log} 2>&1;
        xz -2 {params.uncompressed_alignment};
        xz -2 {params.outdir}/{params.basename}*.fasta
        """

rule diagnostic:
    message: "Scanning metadata {input.metadata} for problematic sequences. Removing sequences with >{params.clock_filter} deviation from the clock and with more than {params.snp_clusters}."
    input:
        metadata = "results/sanitized_metadata_{origin}.tsv.xz"
    output:
        to_exclude = "results/to-exclude_{origin}.txt"
    params:
        clock_filter = 20,
        snp_clusters = 1,
        rare_mutations = 100,
        clock_plus_rare = 100,
    log:
        "logs/diagnostics_{origin}.txt"
    benchmark:
        "benchmarks/diagnostics_{origin}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --metadata {input.metadata} \
            --clock-filter {params.clock_filter} \
            --rare-mutations {params.rare_mutations} \
            --clock-plus-rare {params.clock_plus_rare} \
            --snp-clusters {params.snp_clusters} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

def _collect_exclusion_files(wildcards):
    # This rule creates a per-input exclude file for `rule filter`. This file contains one or both of the following:
    # (1) a config-defined exclude file
    # (2) a dynamically created file (`rule diagnostic`) which scans the alignment for potential errors
    # The second file is optional - it may be opted out via config → skip_diagnostics
    # If the input starting point is "masked" then we also ignore the second file, as the alignment is not available
    if config["filter"].get(wildcards["origin"], {}).get("skip_diagnostics", False):
        return [ config["files"]["exclude"] ]
    if "masked" in config["inputs"][wildcards["origin"]]:
        return [ config["files"]["exclude"] ]
    return [ config["files"]["exclude"], f"results/to-exclude_{wildcards['origin']}.txt" ]

rule mask:
    message:
        """
        Mask bases in alignment {input.alignment}
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = lambda w: _get_path_for_input("aligned", w.origin)
    output:
        alignment = "results/masked_{origin}.fasta.xz"
    log:
        "logs/mask_{origin}.txt"
    benchmark:
        "benchmarks/mask_{origin}.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"],
        mask_sites = config["mask"]["mask_sites"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --mask-terminal-gaps \
            --output /dev/stdout | xz -c -2 > {output.alignment} 2> {log}
        """

rule filter:
    message:
        """
        Filtering alignment {input.sequences} -> {output.sequences}
          - excluding strains in {input.exclude}
          - including strains in {input.include}
          - min length: {params.min_length}
        """
    input:
        sequences = lambda wildcards: _get_path_for_input("masked", wildcards.origin),
        metadata = "results/sanitized_metadata_{origin}.tsv.xz",
        # TODO - currently the include / exclude files are not input (origin) specific, but this is possible if we want
        include = config["files"]["include"],
        exclude = _collect_exclusion_files,
    output:
        sequences = "results/filtered_{origin}.fasta.xz"
    log:
        "logs/filtered_{origin}.txt"
    benchmark:
        "benchmarks/filter_{origin}.txt"
    params:
        min_length = lambda wildcards: _get_filter_value(wildcards, "min_length"),
        exclude_where = lambda wildcards: _get_filter_value(wildcards, "exclude_where"),
        min_date = lambda wildcards: _get_filter_value(wildcards, "min_date"),
        ambiguous = lambda wildcards: f"--exclude-ambiguous-dates-by {_get_filter_value(wildcards, 'exclude_ambiguous_dates_by')}" if _get_filter_value(wildcards, "exclude_ambiguous_dates_by") else "",
        date = (date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d"),
        intermediate_output=lambda wildcards, output: Path(output.sequences).with_suffix("")
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --max-date {params.date} \
            --min-date {params.min_date} \
            {params.ambiguous} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --output {params.intermediate_output} 2>&1 | tee {log};
        xz -2 {params.intermediate_output}
        """

def _get_subsampling_settings(wildcards):
    # Allow users to override default subsampling with their own settings keyed
    # by location type and name. For example, "region_europe" or
    # "country_iceland". Otherwise, default to settings for the location type.
    subsampling_scheme = _get_subsampling_scheme_by_build_name(wildcards.build_name)

    # When there is no well-defined subsampling scheme, default to using all
    # available samples.
    if subsampling_scheme not in config["subsampling"]:
        print(
            f"WARNING: No valid subsampling scheme is defined for build '{wildcards.build_name}'.",
            "Skipping subsampling and using all available samples.",
            file=sys.stderr
        )
        subsampling_scheme = "all"

    subsampling_settings = config["subsampling"][subsampling_scheme]

    if hasattr(wildcards, "subsample"):
        subsampling_settings = subsampling_settings[wildcards.subsample]

        # If users have supplied both `max_sequences` and `seq_per_group`, we
        # throw an error instead of assuming the user prefers one setting over
        # another by default.
        if subsampling_settings.get("max_sequences") and subsampling_settings.get("seq_per_group"):
            raise Exception(f"The subsampling scheme '{subsampling_scheme}' for build '{wildcards.build_name}' defines both `max_sequences` and `seq_per_group`, but these arguments are mutually exclusive. If you didn't define both of these settings, this conflict could be caused by using the same subsampling scheme name as a default scheme. In this case, rename your subsampling scheme, '{subsampling_scheme}', to a unique name (e.g., 'custom_{subsampling_scheme}') and run the workflow again.")

        # If users have defined `group_by` but supplied neither `max_sequences` nor `seq_per_group`, we
        # throw an error because the subsampling rule will still group by one or
        # more fields and the lack of limits on this grouping could produce
        # unexpected behavior.
        if subsampling_settings.get("group_by") and not subsampling_settings.get("max_sequences") and not subsampling_settings.get("seq_per_group"):
            raise Exception(f"The subsampling scheme '{subsampling_scheme}' for build '{wildcards.build_name}' must define `max_sequences` or `seq_per_group`.")

    return subsampling_settings


rule combine_sequences_for_subsampling:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling.
        """
    input:
        lambda w: [_get_path_for_input("filtered", origin) for origin in config.get("inputs", {})]
    output:
        "results/combined_sequences_for_subsampling.fasta.xz"
    benchmark:
        "benchmarks/combine_sequences_for_subsampling.txt"
    conda: config["conda_environment"]
    params:
        error_on_duplicate_strains="--error-on-duplicate-strains" if not config.get("combine_sequences_for_subsampling", {}).get("warn_about_duplicates") else "",
        strain_prefixes=config["strip_strain_prefixes"],
    shell:
        """
        python3 scripts/sanitize_sequences.py \
                --sequences {input} \
                --strip-prefixes {params.strain_prefixes:q} \
                {params.error_on_duplicate_strains} \
                --output /dev/stdout \
                | xz -c -2 > {output}
        """

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = _get_unified_alignment
    output:
        sequence_index = "results/combined_sequence_index.tsv.xz"
    log:
        "logs/index_sequences.txt"
    benchmark:
        "benchmarks/index_sequences.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index} 2>&1 | tee {log}
        """

rule extract_subsampling_scheme:
    message:
        """
        Extracting subsampling scheme for build "{wildcards.build_name}" into its own YAML file
        """
    output:
        scheme="results/{build_name}/subsampling_scheme.yaml",
    run:
        # Note that the syntax is slightly different between the YAML which `./scripts/subsample.py`
        # expects (this script is a precursor for `augur subsample`) and the way in which we write
        # subsampling definitions in our `builds.yaml` file.
        # (1) wildcards must be filled in
        # (2) we use `augur filter`-like syntax, e.g. "sequences-per-group" not "seq_per_group"
        import yaml
        import json
        import shlex
        scheme = _get_subsampling_settings(wildcards)

        # fill in templates, e.g `{country}`, with build-specific data
        # (prior art existed in the function _get_specific_subsampling_setting)
        build = config["builds"][wildcards.build_name]
        for sample_description in scheme.values():
            for key,value in sample_description.items():
                if isinstance(value, str):
                    sample_description[key]=value.format(**build)

        # we allow "no_subsampling" via a boolean in _each_ sample definition
        # and pass each sample through the subsampling process, as we used to
        # (this should be optimised via snakemake logic)
        for sample_name, sample_dict in scheme.items():
            if sample_dict.get("no_subsampling", False) == True:
                scheme[sample_name] = {}

        # map ncov-specific syntax into general subsampling syntax
        # (prior art existed in the function _get_specific_subsampling_setting)
        for sample in scheme.values():
            if "group_by" in sample:
                sample["group-by"] = shlex.split(sample['group_by'])
                del sample["group_by"]
            if "seq_per_group" in sample:
                sample["sequences-per-group"] = sample['seq_per_group']
                del sample["seq_per_group"]
            if "max_sequences" in sample:
                sample["subsample-max-sequences"] = sample['max_sequences']
                del sample["max_sequences"]
            if "exclude_ambiguous_dates_by" in sample:
                sample["exclude-ambiguous-dates-by"] = sample["exclude_ambiguous_dates_by"]
                del sample["exclude_ambiguous_dates_by"]

            # Certain settings required that the argument was defined in the value,
            # e.g. max_date: "--max-date 2020-01-02", which we remove here
            if "query" in sample: # ncov syntax included the `--query` argument
                sample["query"] = sample['query'].lstrip("--query ").strip('"').strip("'")
            if "min_date" in sample: # ncov syntax included `--min-date` in the value
                sample["min-date"] = sample["min_date"].lstrip("--min-date ")
                del sample["min_date"]
            if "max_date" in sample: # ncov syntax included `--max-date` in the value
                sample["max-date"] = sample["max_date"].lstrip("--max-date ")
                del sample["max_date"]

            # the include/exclude ncov subsampling settings meant {include,exclude}-where
            # and also included the argument in their value, which we remove here
            if "exclude" in sample:
                sample["exclude-where"] = shlex.split(sample['exclude'])[1:]
                del sample["exclude"]
            if "include" in sample:
                sample["include-where"] = shlex.split(sample['include'])[1:]
                del sample["include"]

            if "sampling_scheme" in sample: # ncov syntax for expressing additional filter arguments
                if sample["sampling_scheme"]=="--probabilistic-sampling":
                    sample["probabilistic-sampling"] = True
                elif sample["sampling_scheme"]=="--no-probabilistic-sampling":
                    sample["no-probabilistic-sampling"] = True
                del sample["sampling_scheme"]

        with open(output.scheme, 'w', encoding='utf-8') as fh:
            # to avoid a YAML with type annotations, we roundtrip through a JSON
            yaml.dump(json.loads(json.dumps(scheme)), fh, default_flow_style=False)

rule subsample:
    message:
        """
        Subsampling using our subsampling script.
        """
    input:
        scheme="results/{build_name}/subsampling_scheme.yaml",
        metadata=_get_unified_metadata,
        alignment=_get_unified_alignment,
        alignment_index=rules.index_sequences.output.sequence_index,
        reference = config["files"]["alignment_reference"],
        include = config["files"]["include"],
        exclude = config["files"]["exclude"]
    output:
        intermediates=directory("results/{build_name}/subsamples/"),
        sequences = "results/{build_name}/{build_name}_subsampled_sequences.fasta.xz",
        metadata = "results/{build_name}/{build_name}_subsampled_metadata.tsv.xz",
        # output_log = "results/{build_name}/subsamples/subsampled.tsv",
    log:
        "logs/subsample_{build_name}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}.txt"
    conda: config["conda_environment"]
    threads: 4
    resources:
        mem_mb=12000
    shell:
        """
        python3 scripts/subsample.py \
            --scheme {input.scheme} \
            --metadata {input.metadata} \
            --alignment {input.alignment} \
            --alignment-index {input.alignment_index} \
            --reference {input.reference} \
            --include-strains-file {input.include} \
            --exclude-strains-file {input.exclude} \
            --output-dir {output.intermediates} \
            --output-fasta {output.sequences} \
            --output-metadata {output.metadata} 2>&1 | tee {log}
        """

rule build_align:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = rules.subsample.output.sequences,
        genemap = config["files"]["annotation"],
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "results/{build_name}/aligned.fasta",
        insertions = "results/{build_name}/insertions.tsv",
        translations = expand("results/{{build_name}}/translations/aligned.gene.{gene}.fasta", gene=config.get('genes', ['S']))
    params:
        outdir = "results/{build_name}/translations",
        genes = ','.join(config.get('genes', ['S'])),
        basename = "aligned"
    log:
        "logs/align_{build_name}.txt"
    benchmark:
        "benchmarks/align_{build_name}.txt"
    conda: config["conda_environment"]
    threads: 8
    resources:
        mem_mb=3000
    shell:
        """
        xz -c -d {input.sequences} | nextalign \
            --jobs={threads} \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --genes {params.genes} \
            --sequences /dev/stdin \
            --output-dir {params.outdir} \
            --output-basename {params.basename} \
            --output-fasta {output.alignment} \
            --output-insertions {output.insertions} > {log} 2>&1
        """

rule compress_build_align:
    message:
        """Compressing {input.alignment}"""
    input:
        alignment = "results/{build_name}/aligned.fasta"
    output:
        alignment = "results/{build_name}/aligned.fasta.xz"
    benchmark:
        "benchmarks/compress_build_align_{build_name}.txt"
    conda:
        config["conda_environment"]
    log:
        "logs/compress_build_align_{build_name}.txt"
    shell:
        """
        xz -c {input} > {output} 2> {log}
        """

if "run_pangolin" in config and config["run_pangolin"]:
    rule run_pangolin:
        message:
            """
            Running pangolin to assign lineage labels to samples. Includes putative lineage definitions by default.
            Please remember to update your installation of pangolin regularly to ensure the most up-to-date classifications.
            """
        input:
            alignment = rules.build_align.output.alignment,
        output:
            lineages = "results/{build_name}/pangolineages.csv",
        params:
            outdir = "results/{build_name}",
            csv_outfile = "pangolineages.csv",
            node_data_outfile = "pangolineages.json"
        log:
            "logs/pangolin_{build_name}.txt"
        conda: config["conda_environment"]
        threads: 1
        resources:
            mem_mb=3000
        benchmark:
            "benchmarks/pangolineages_{build_name}.txt"
        shell: ## once pangolin fully supports threads, add `--threads {threads}` to the below (existing pango cli param)
            """
            pangolin {input.alignment}\
                --outdir {params.outdir} \
                --outfile {params.csv_outfile} 2>&1 | tee {log}\
            """

    rule make_pangolin_node_data:
        input:
            lineages = rules.run_pangolin.output.lineages
        output:
            node_data = "results/{build_name}/pangolineages.json"
        log:
            "logs/pangolin_export_{build_name}.txt"
        conda: config["conda_environment"]
        resources:
            mem_mb=3000
        benchmark:
            "benchmarks/make_pangolin_node_data_{build_name}.txt"
        shell:
            """
            python3 scripts/make_pangolin_node_data.py \
            --pangolineages {input.lineages} \
            --node_data_outfile {output.node_data} 2>&1 | tee {log}\
            """

# TODO: This will probably not work for build names like "country_usa" where we need to know the country is "USA".
rule adjust_metadata_regions:
    message:
        """
        Adjusting metadata for build '{wildcards.build_name}'
        """
    input:
        metadata="results/{build_name}/{build_name}_subsampled_metadata.tsv.xz",
    output:
        metadata = "results/{build_name}/metadata_adjusted.tsv.xz"
    params:
        # Default to a "global" region if none is defined. The adjust metadata
        # script will not modify the metadata if the region is "global".
        region = lambda wildcards: config["builds"][wildcards.build_name].get("region", "global")
    log:
        "logs/adjust_metadata_regions_{build_name}.txt"
    benchmark:
        "benchmarks/adjust_metadata_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/adjust_regional_meta.py \
            --region {params.region:q} \
            --metadata {input.metadata} \
            --output {output.metadata} 2>&1 | tee {log}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.build_align.output.alignment
    output:
        tree = "results/{build_name}/tree_raw.nwk"
    params:
        args = lambda w: config["tree"].get("tree-builder-args","") if "tree" in config else ""
    log:
        "logs/tree_{build_name}.txt"
    benchmark:
        "benchmarks/tree_{build_name}.txt"
    threads: 8
    resources:
        # Multiple sequence alignments can use up to 40 times their disk size in
        # memory, especially for larger alignments.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 40 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.build_align.output.alignment,
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        tree = "results/{build_name}/tree.nwk",
        node_data = "results/{build_name}/branch_lengths.json"
    log:
        "logs/refine_{build_name}.txt"
    benchmark:
        "benchmarks/refine_{build_name}.txt"
    threads: 1
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    params:
        root = config["refine"]["root"],
        clock_rate = config["refine"]["clock_rate"],
        clock_std_dev = config["refine"]["clock_std_dev"],
        coalescent = config["refine"]["coalescent"],
        date_inference = config["refine"]["date_inference"],
        divergence_unit = config["refine"]["divergence_unit"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"],
        keep_polytomies = "--keep-polytomies" if config["refine"].get("keep_polytomies", False) else "",
        timetree = "" if config["refine"].get("no_timetree", False) else "--timetree"
    conda: config["conda_environment"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            {params.timetree} \
            {params.keep_polytomies} \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd} 2>&1 | tee {log}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.build_align.output.alignment
    output:
        node_data = "results/{build_name}/nt_muts.json"
    log:
        "logs/ancestral_{build_name}.txt"
    benchmark:
        "benchmarks/ancestral_{build_name}.txt"
    params:
        inference = config["ancestral"]["inference"]
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous 2>&1 | tee {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = config["files"]["reference"]
    output:
        node_data = "results/{build_name}/aa_muts.json"
    log:
        "logs/translate_{build_name}.txt"
    benchmark:
        "benchmarks/translate_{build_name}.txt"
    resources:
        # Memory use scales primarily with size of the node data.
        mem_mb=lambda wildcards, input: 3 * int(input.node_data.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule aa_muts_explicit:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        translations = lambda w: rules.build_align.output.translations
    output:
        node_data = "results/{build_name}/aa_muts_explicit.json",
        translations = expand("results/{{build_name}}/translations/aligned.gene.{gene}_withInternalNodes.fasta", gene=config.get('genes', ['S']))
    params:
        genes = config.get('genes', 'S')
    log:
        "logs/aamuts_{build_name}.txt"
    benchmark:
        "benchmarks/aamuts_{build_name}.txt"
    resources:
        # Multiple sequence alignments can use up to 15 times their disk size in
        # memory.
        # Note that Snakemake >5.10.0 supports input.size_mb to avoid converting from bytes to MB.
        mem_mb=lambda wildcards, input: 15 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/explicit_translation.py \
            --tree {input.tree} \
            --translations {input.translations:q} \
            --genes {params.genes} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule build_mutation_summary:
    message: "Summarizing {input.alignment}"
    input:
        alignment = rules.build_align.output.alignment,
        insertions = rules.build_align.output.insertions,
        translations = rules.build_align.output.translations,
        reference = config["files"]["alignment_reference"],
        genemap = config["files"]["annotation"]
    output:
        mutation_summary = "results/{build_name}/mutation_summary.tsv"
    log:
        "logs/mutation_summary_{build_name}.txt"
    params:
        outdir = "results/{build_name}/translations",
        basename = "aligned"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mutation_summary.py \
            --alignment {input.alignment} \
            --insertions {input.insertions} \
            --directory {params.outdir} \
            --basename {params.basename} \
            --reference {input.reference} \
            --genemap {input.genemap} \
            --output {output.mutation_summary} 2>&1 | tee {log}
        """

rule distances:
    input:
        tree = rules.refine.output.tree,
        alignments = "results/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
        distance_maps = ["defaults/distance_maps/S1.json"]
    params:
        genes = 'S',
        comparisons = ['root'],
        attribute_names = ['S1_mutations']
    output:
        node_data = "results/{build_name}/distances.json"
    conda:
        config["conda_environment"]
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output}
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        node_data = "results/{build_name}/traits.json"
    log:
        "logs/traits_{build_name}.txt"
    benchmark:
        "benchmarks/traits_{build_name}.txt"
    params:
        columns = _get_trait_columns_by_wildcards,
        sampling_bias_correction = _get_sampling_bias_correction_for_wildcards
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} 2>&1 | tee {log}
        """

def _get_clade_files(wildcards):
    if "subclades" in config["builds"][wildcards.build_name]:
        return [config["files"]["clades"], config["builds"][wildcards.build_name]["subclades"]]
    else:
        return config["files"]["clades"]

rule clade_files:
    input:
        clade_files = _get_clade_files
    output:
        "results/{build_name}/clades.tsv"
    benchmark:
        "benchmarks/clade_files_{build_name}.txt"
    shell:
        '''
        cat {input.clade_files} > {output}
        '''

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = rules.clade_files.output
    output:
        clade_data = "results/{build_name}/clades.json"
    log:
        "logs/clades_{build_name}.txt"
    benchmark:
        "benchmarks/clades_{build_name}.txt"
    resources:
        # Memory use scales primarily with size of the node data.
        mem_mb=lambda wildcards, input: 3 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data} 2>&1 | tee {log}
        """

rule emerging_lineages:
    message: "Adding emerging clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        emerging_lineages = config["files"]["emerging_lineages"],
        clades = config["files"]["clades"]
    output:
        clade_data = "results/{build_name}/temp_emerging_lineages.json"
    log:
        "logs/emerging_lineages_{build_name}.txt"
    benchmark:
        "benchmarks/emerging_lineages_{build_name}.txt"
    resources:
        # Memory use scales primarily with size of the node data.
        mem_mb=lambda wildcards, input: 3 * int(input.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.emerging_lineages} \
            --output-node-data {output.clade_data} 2>&1 | tee {log}
        """

rule rename_emerging_lineages:
    input:
        node_data = rules.emerging_lineages.output.clade_data
    output:
        clade_data = "results/{build_name}/emerging_lineages.json"
    benchmark:
        "benchmarks/rename_emerging_lineages_{build_name}.txt"
    run:
        import json
        with open(input.node_data, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            new_data = {}
            for k,v in d['nodes'].items():
                if "clade_membership" in v:
                    new_data[k] = {"emerging_lineage": v["clade_membership"]}
        with open(output.clade_data, "w") as fh:
            json.dump({"nodes": new_data}, fh, indent=2)


rule colors:
    message: "Constructing colors file"
    input:
        ordering = config["files"]["ordering"],
        color_schemes = config["files"]["color_schemes"],
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        colors = "results/{build_name}/colors.tsv"
    log:
        "logs/colors_{build_name}.txt"
    benchmark:
        "benchmarks/colors_{build_name}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        # Compared to other rules, this rule loads metadata as a pandas
        # DataFrame instead of a dictionary, so it uses much less memory.
        mem_mb=lambda wildcards, input: 5 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1 | tee {log}
        """

rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        node_data = "results/{build_name}/recency.json"
    log:
        "logs/recency_{build_name}.txt"
    benchmark:
        "benchmarks/recency_{build_name}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output} 2>&1 | tee {log}
        """

rule tip_frequencies:
    message: "Estimating censored KDE frequencies for tips"
    input:
        tree = rules.refine.output.tree,
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        tip_frequencies_json = "results/{build_name}/tip-frequencies.json"
    log:
        "logs/tip_frequencies_{build_name}.txt"
    benchmark:
        "benchmarks/tip_frequencies_{build_name}.txt"
    params:
        min_date = _get_min_date_for_frequencies,
        max_date = _get_max_date_for_frequencies,
        pivot_interval = config["frequencies"]["pivot_interval"],
        pivot_interval_units = config["frequencies"]["pivot_interval_units"],
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"]
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_interval} \
            --pivot-interval-units {params.pivot_interval_units} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --output {output.tip_frequencies_json} 2>&1 | tee {log}
        """

rule logistic_growth:
    input:
        tree="results/{build_name}/tree.nwk",
        frequencies="results/{build_name}/tip-frequencies.json",
    output:
        node_data="results/{build_name}/logistic_growth.json"
    benchmark:
        "benchmarks/logistic_growth_{build_name}.txt"
    conda:
        config["conda_environment"]
    log:
        "logs/logistic_growth_{build_name}.txt"
    params:
        method="logistic",
        attribute_name = "logistic_growth",
        delta_pivots=config["logistic_growth"]["delta_pivots"],
        min_tips=config["logistic_growth"]["min_tips"],
        min_frequency=config["logistic_growth"]["min_frequency"],
        max_frequency=config["logistic_growth"]["max_frequency"],
    resources:
        mem_mb=256
    shell:
        """
        python3 scripts/calculate_delta_frequency.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --method {params.method} \
            --delta-pivots {params.delta_pivots} \
            --min-tips {params.min_tips} \
            --min-frequency {params.min_frequency} \
            --max-frequency {params.max_frequency} \
            --attribute-name {params.attribute_name} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule mutational_fitness:
    input:
        tree = "results/{build_name}/tree.nwk",
        alignments = lambda w: rules.aa_muts_explicit.output.translations,
        distance_map = config["files"]["mutational_fitness_distance_map"]
    output:
        node_data = "results/{build_name}/mutational_fitness.json"
    benchmark:
        "benchmarks/mutational_fitness_{build_name}.txt"
    conda:
        config["conda_environment"]
    log:
        "logs/mutational_fitness_{build_name}.txt"
    params:
        genes = ' '.join(config.get('genes', ['S'])),
        compare_to = "root",
        attribute_name = "mutational_fitness"
    conda:
        config["conda_environment"],
    resources:
        mem_mb=2000
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.compare_to} \
            --attribute-name {params.attribute_name} \
            --map {input.distance_map} \
            --output {output} 2>&1 | tee {log}
        """

rule calculate_epiweeks:
    input:
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        node_data="results/{build_name}/epiweeks.json",
    benchmark:
        "benchmarks/calculate_epiweeks_{build_name}.txt",
    conda:
        config["conda_environment"],
    log:
        "logs/calculate_epiweeks_{build_name}.txt",
    shell:
        """
        python3 scripts/calculate_epiweek.py \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

def export_title(wildcards):
    # TODO: maybe we could replace this with a config entry for full/human-readable build name?
    location_name = wildcards.build_name

    # If specified in config file generally, or in a config file build
    if "title" in config["builds"][location_name]:
        return config["builds"][location_name]["title"]
    elif "title" in config:
        return config["title"]

    # Else return an auto-generated title
    if not location_name:
        return "Genomic epidemiology of novel coronavirus"
    elif location_name == "global":
        return "Genomic epidemiology of novel coronavirus - Global subsampling"
    else:
        location_title = location_name.replace("-", " ").title()
        return f"Genomic epidemiology of novel coronavirus - {location_title}-focused subsampling"

def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.rename_emerging_lineages.output.clade_data,
        rules.clades.output.clade_data,
        rules.recency.output.node_data,
        rules.traits.output.node_data,
        rules.logistic_growth.output.node_data,
        rules.mutational_fitness.output.node_data,
        rules.aa_muts_explicit.output.node_data,
        rules.distances.output.node_data,
        rules.calculate_epiweeks.output.node_data
    ]

    if "run_pangolin" in config and config["run_pangolin"]:
        inputs.append(rules.make_pangolin_node_data.output.node_data)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]

    return inputs

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
        node_data = _get_node_data_by_wildcards,
        auspice_config = lambda w: config["builds"][w.build_name]["auspice_config"] if "auspice_config" in config["builds"][w.build_name] else config["files"]["auspice_config"],
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"][w.build_name] else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"],
        description = lambda w: config["builds"][w.build_name]["description"] if "description" in config["builds"][w.build_name] else config["files"]["description"]
    output:
        auspice_json = "results/{build_name}/ncov_with_accessions.json",
        root_sequence_json = "results/{build_name}/ncov_with_accessions_root-sequence.json"
    log:
        "logs/export_{build_name}.txt"
    benchmark:
        "benchmarks/export_{build_name}.txt"
    params:
        title = export_title
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=12000
    conda: config["conda_environment"]
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule add_branch_labels:
    message: "Adding custom branch labels to the Auspice JSON"
    input:
        auspice_json = rules.export.output.auspice_json,
        emerging_clades = rules.emerging_lineages.output.clade_data
    output:
        auspice_json = "results/{build_name}/ncov_with_branch_labels.json"
    log:
        "logs/add_branch_labels{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 ./scripts/add_branch_labels.py \
            --input {input.auspice_json} \
            --emerging-clades {input.emerging_clades} \
            --output {output.auspice_json}
        """

rule include_hcov19_prefix:
    message: "Rename strains to include hCoV-19/ prefix"
    input:
        auspice_json = rules.add_branch_labels.output.auspice_json,
        tip_frequencies = rules.tip_frequencies.output.tip_frequencies_json
    output:
        auspice_json = "results/{build_name}/ncov_with_hcov19_prefix.json",
        tip_frequencies = "results/{build_name}/tip-frequencies_with_hcov19_prefix.json"
    log:
        "logs/include_hcov19_prefix{build_name}.txt"
    conda: config["conda_environment"]
    params:
        prefix = lambda w: "hCoV-19/" if config.get("include_hcov19_prefix", False) else ""
    shell:
        """
        python3 ./scripts/include_prefix.py \
            --input-auspice {input.auspice_json} \
            --input-tip-frequencies {input.tip_frequencies} \
            --prefix {params.prefix} \
            --output-auspice {output.auspice_json} \
            --output-tip-frequencies {output.tip_frequencies}
        """

rule finalize:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = lambda w: rules.include_hcov19_prefix.output.auspice_json,
        frequencies = rules.include_hcov19_prefix.output.tip_frequencies,
        root_sequence_json = rules.export.output.root_sequence_json
    output:
        auspice_json = f"auspice/{config['auspice_json_prefix']}_{{build_name}}.json",
        tip_frequency_json = f"auspice/{config['auspice_json_prefix']}_{{build_name}}_tip-frequencies.json",
        root_sequence_json = f"auspice/{config['auspice_json_prefix']}_{{build_name}}_root-sequence.json"
    log:
        "logs/fix_colorings_{build_name}.txt"
    benchmark:
        "benchmarks/fix_colorings_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log} &&
        cp {input.frequencies} {output.tip_frequency_json} &&
        cp {input.root_sequence_json} {output.root_sequence_json}
        """
