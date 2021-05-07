rule sanitize_metadata:
    input:
        metadata=lambda wildcards: _get_path_for_input("metadata", wildcards.origin)
    output:
        metadata="pre-processed/sanitized_metadata{origin}.tsv"
    benchmark:
        "benchmarks/sanitize_metadata{origin}.txt"
    conda:
        config["conda_environment"]
    log:
        "logs/sanitize_metadata{origin}.txt"
    params:
        strain_prefixes=config["strip_strain_prefixes"],
    shell:
        """
        python3 scripts/sanitize_metadata.py \
            --metadata {input.metadata} \
            --strip-prefixes {params.strain_prefixes:q} \
            --output {output.metadata}
        """


rule combine_input_metadata:
    # this rule is intended to be run _only_ if we have defined multiple inputs ("origins")
    message:
        """
        Combining metadata files {input.metadata} -> {output.metadata} and adding columns to represent origin
        """
    input:
        metadata=expand("pre-processed/sanitized_metadata{origin}.tsv", origin=[f"_{origin}" if origin else "" for origin in config.get("inputs", "")]),
    output:
        metadata = "pre-processed/combined_metadata.tsv"
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

if "use_nextalign" in config and config["use_nextalign"]:
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
            alignment = "pre-processed/aligned{origin}.fasta",
            insertions = "pre-processed/insertions{origin}.tsv",
            translations = expand("pre-processed/translations/seqs{{origin}}.gene.{gene}.fasta", gene=config.get('genes', ['S']))
        params:
            outdir = "pre-processed/translations",
            genes = ','.join(config.get('genes', ['S'])),
            basename = "seqs{origin}",
            strain_prefixes=config["strip_strain_prefixes"],
        log:
            "logs/align{origin}.txt"
        benchmark:
            "benchmarks/align{origin}.txt"
        conda: config["conda_environment"]
        threads: 8
        resources:
            mem_mb=3000
        shell:
            """
            python3 scripts/sanitize_sequences.py \
                --sequences {input.sequences} \
                --strip-prefixes {params.strain_prefixes:q} \
                --output /dev/stdout \
                | nextalign \
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
else:
    rule align:
        message:
            """
            Aligning sequences from {input.sequences} to {input.reference}
            - gaps relative to reference are considered real
            """
        input:
            sequences = lambda wildcards: _get_path_for_input("sequences", wildcards.origin),
            reference = config["files"]["alignment_reference"]
        output:
            alignment = "pre-processed/aligned{origin}.fasta"
        log:
            "logs/align{origin}.txt"
        benchmark:
            "benchmarks/align{origin}.txt"
        threads: 16
        conda: config["conda_environment"]
        params:
            strain_prefixes=config["strip_strain_prefixes"],
        shell:
            """
            python3 scripts/sanitize_sequences.py \
                --sequences {input.sequences} \
                --strip-prefixes {params.strain_prefixes:q} \
                --output /dev/stdout \
                | mafft \
                --auto \
                --thread {threads} \
                --keeplength \
                --addfragments \
                /dev/stdin \
                {input.reference} > {output} 2> {log}
            """

rule diagnostic:
    message: "Scanning aligned sequences {input.alignment} for problematic sequences"
    input:
        alignment = lambda wildcards: _get_path_for_input("aligned", wildcards.origin),
        metadata = "pre-processed/sanitized_metadata{origin}.tsv",
        reference = config["files"]["reference"]
    output:
        diagnostics = "pre-processed/sequence-diagnostics{origin}.tsv",
        flagged = "pre-processed/flagged-sequences{origin}.tsv",
        to_exclude = "pre-processed/to-exclude{origin}.txt"
    log:
        "logs/diagnostics{origin}.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"]
    benchmark:
        "benchmarks/diagnostics{origin}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --output-flagged {output.flagged} \
            --output-diagnostics {output.diagnostics} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

def _collect_exclusion_files(wildcards):
    # Note that we _always_ exclude the sequences from the (config-defined) exclude file
    # As well as the sequences flagged by the diagnostic step.
    # Note that we can skip the diagnostic step on a per-input (per-origin) basis.
    exclude_files = [ config["files"]["exclude"] ]
    if not config["filter"].get(_trim_origin(wildcards["origin"]), {}).get("skip_diagnostics", False):
        exclude_files.append(_get_path_for_input("to-exclude", wildcards.origin))
    return exclude_files

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
        alignment = "pre-processed/masked{origin}.fasta"
    log:
        "logs/mask{origin}.txt"
    benchmark:
        "benchmarks/mask{origin}.txt"
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
            --output {output.alignment} 2>&1 | tee {log}
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
        metadata = "pre-processed/sanitized_metadata{origin}.tsv",
        # TODO - currently the include / exclude files are not input (origin) specific, but this is possible if we want
        include = config["files"]["include"],
        exclude = _collect_exclusion_files,
    output:
        sequences = "pre-processed/filtered{origin}.fasta"
    log:
        "logs/filtered{origin}.txt"
    benchmark:
        "benchmarks/filter{origin}.txt"
    params:
        min_length = lambda wildcards: _get_filter_value(wildcards, "min_length"),
        exclude_where = lambda wildcards: _get_filter_value(wildcards, "exclude_where"),
        min_date = lambda wildcards: _get_filter_value(wildcards, "min_date"),
        ambiguous = lambda wildcards: f"--exclude-ambiguous-dates-by {_get_filter_value(wildcards, 'exclude_ambiguous_dates_by')}" if _get_filter_value(wildcards, "exclude_ambiguous_dates_by") else "",
        date = (date.today() + datetime.timedelta(days=1)).strftime("%Y-%m-%d")
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
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
            --output {output.sequences} 2>&1 | tee {log}
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


def get_priorities(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return f"builds/{wildcards.build_name}/priorities_{subsampling_settings['priorities']['focus']}.tsv"
    else:
        # TODO: find a way to make the list of input files depend on config
        return config["files"]["include"]


def get_priority_argument(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return "--priority " + get_priorities(wildcards)
    else:
        return ""


def _get_specific_subsampling_setting(setting, optional=False):
    # Note -- this function contains a lot of conditional logic because
    # we have the situation where some config options must define the
    # augur argument in their value, and some must not. For instance:
    # subsamplingScheme -> sampleName -> group_by: year                            (`--group-by` is _not_ part of this value)
    #                                 -> exclude: "--exclude-where 'country=USA'"  (`--exclude-where` IS part of this value)
    # Since there are a lot of subsampling schemes out there, backwards compatability
    # is important!                                 james hadfield, feb 2021
    def _get_setting(wildcards):
        if optional:
            value = _get_subsampling_settings(wildcards).get(setting, "")
        else:
            value = _get_subsampling_settings(wildcards)[setting]

        if isinstance(value, str):
            # Load build attributes including geographic details about the
            # build's region, country, division, etc. as needed for subsampling.
            build = config["builds"][wildcards.build_name]
            value = value.format(**build)
            if value !="":
                if setting == 'exclude_ambiguous_dates_by':
                    value = f"--exclude-ambiguous-dates-by {value}"
                elif setting == 'group_by':
                    value = f"--group-by {value}"
        elif value is not None:
            # If is 'seq_per_group' or 'max_sequences' build subsampling setting,
            # need to return the 'argument' for augur
            if setting == 'seq_per_group':
                value = f"--sequences-per-group {value}"
            elif setting == 'max_sequences':
                value = f"--subsample-max-sequences {value}"

            return value
        else:
            value = ""

        # Check format strings that haven't been resolved.
        if re.search(r'\{.+\}', value):
            raise Exception(f"The parameters for the subsampling scheme '{wildcards.subsample}' of build '{wildcards.build_name}' reference build attributes that are not defined in the configuration file: '{value}'. Add these build attributes to the appropriate configuration file and try again.")

        return value

    return _get_setting


rule combine_sequences_for_subsampling:
    # Similar to rule combine_input_metadata, this rule should only be run if multiple inputs are being used (i.e. multiple origins)
    message:
        """
        Combine and deduplicate aligned & filtered FASTAs from multiple origins in preparation for subsampling.
        """
    input:
        lambda w: [_get_path_for_input("filtered", f"_{origin}") for origin in config.get("inputs", {})]
    output:
        "pre-processed/combined_sequences_for_subsampling.fasta"
    benchmark:
        "benchmarks/combine_sequences_for_subsampling.txt"
    conda: config["conda_environment"]
    params:
        warn_about_duplicates="--warn-about-duplicates" if config.get("combine_sequences_for_subsampling", {}).get("warn_about_duplicates") else "",
        strain_prefixes=config["strip_strain_prefixes"],
    shell:
        """
        python3 scripts/sanitize_sequences.py \
                --sequences {input} \
                --strip-prefixes {params.strain_prefixes:q} \
                --output /dev/stdout \
                | python3 scripts/combine-and-dedup-fastas.py \
            --input /dev/stdin \
            {params.warn_about_duplicates} \
            --output {output}
        """

rule index_sequences:
    message:
        """
        Index sequence composition for faster filtering.
        """
    input:
        sequences = _get_unified_alignment
    output:
        sequence_index = "pre-processed/combined_sequence_index.tsv"
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

rule subsample:
    message:
        """
        Subsample all sequences by '{wildcards.subsample}' scheme for build '{wildcards.build_name}' with the following parameters:

         - group by: {params.group_by}
         - sequences per group: {params.sequences_per_group}
         - subsample max sequences: {params.subsample_max_sequences}
         - min-date: {params.min_date}
         - max-date: {params.max_date}
         - {params.exclude_ambiguous_dates_argument}
         - exclude: {params.exclude_argument}
         - include: {params.include_argument}
         - query: {params.query_argument}
         - priority: {params.priority_argument}
        """
    input:
        sequences = _get_unified_alignment,
        metadata = _get_unified_metadata,
        sequence_index = rules.index_sequences.output.sequence_index,
        include = config["files"]["include"],
        priorities = get_priorities,
        exclude = config["files"]["exclude"]
    output:
        sequences = "builds/{build_name}/sample-{subsample}.fasta",
        strains="builds/{build_name}/sample-{subsample}.txt",
    log:
        "logs/subsample_{build_name}_{subsample}.txt"
    benchmark:
        "benchmarks/subsample_{build_name}_{subsample}.txt"
    params:
        group_by = _get_specific_subsampling_setting("group_by", optional=True),
        sequences_per_group = _get_specific_subsampling_setting("seq_per_group", optional=True),
        subsample_max_sequences = _get_specific_subsampling_setting("max_sequences", optional=True),
        sampling_scheme = _get_specific_subsampling_setting("sampling_scheme", optional=True),
        exclude_argument = _get_specific_subsampling_setting("exclude", optional=True),
        include_argument = _get_specific_subsampling_setting("include", optional=True),
        query_argument = _get_specific_subsampling_setting("query", optional=True),
        exclude_ambiguous_dates_argument = _get_specific_subsampling_setting("exclude_ambiguous_dates_by", optional=True),
        min_date = _get_specific_subsampling_setting("min_date", optional=True),
        max_date = _get_specific_subsampling_setting("max_date", optional=True),
        priority_argument = get_priority_argument
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --include {input.include} \
            --exclude {input.exclude} \
            {params.min_date} \
            {params.max_date} \
            {params.exclude_argument} \
            {params.include_argument} \
            {params.query_argument} \
            {params.exclude_ambiguous_dates_argument} \
            {params.priority_argument} \
            {params.group_by} \
            {params.sequences_per_group} \
            {params.subsample_max_sequences} \
            {params.sampling_scheme} \
            --output {output.sequences} \
            --output-strains {output.strains} 2>&1 | tee {log}
        """

rule proximity_score:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for build '{wildcards.build_name}'.
        """
    input:
        alignment = _get_unified_alignment,
        reference = config["files"]["alignment_reference"],
        focal_alignment = "builds/{build_name}/sample-{focus}.fasta"
    output:
        proximities = "builds/{build_name}/proximity_{focus}.tsv"
    log:
        "logs/subsampling_proximity_{build_name}_{focus}.txt"
    benchmark:
        "benchmarks/proximity_score_{build_name}_{focus}.txt"
    params:
        chunk_size=10000,
        ignore_seqs = config['refine']['root']
    resources:
        # Memory scales at ~0.15 MB * chunk_size (e.g., 0.15 MB * 10000 = 1.5GB).
        mem_mb=4000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/get_distance_to_focal_set.py \
            --reference {input.reference} \
            --alignment {input.alignment} \
            --focal-alignment {input.focal_alignment} \
            --ignore-seqs {params.ignore_seqs} \
            --chunk-size {params.chunk_size} \
            --output {output.proximities} 2>&1 | tee {log}
        """

rule priority_score:
    input:
        proximity = rules.proximity_score.output.proximities,
        sequence_index = rules.index_sequences.output.sequence_index,
    output:
        priorities = "builds/{build_name}/priorities_{focus}.tsv"
    benchmark:
        "benchmarks/priority_score_{build_name}_{focus}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/priorities.py \
            --sequence-index {input.sequence_index} \
            --proximities {input.proximity} \
            --output {output.priorities} 2>&1 | tee {log}
        """


def _get_subsampled_files(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    return [
        f"builds/{wildcards.build_name}/sample-{subsample}.txt"
        for subsample in subsampling_settings
    ]

rule combine_samples:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        sequences=_get_unified_alignment,
        sequence_index=rules.index_sequences.output.sequence_index,
        metadata=_get_unified_metadata,
        include=_get_subsampled_files,
    output:
        sequences = "builds/{build_name}/{build_name}_subsampled_sequences.fasta",
        metadata = "builds/{build_name}/{build_name}_subsampled_metadata.tsv"
    log:
        "logs/subsample_regions_{build_name}.txt"
    benchmark:
        "benchmarks/subsample_regions_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.include} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} 2>&1 | tee {log}
        """

if "use_nextalign" in config and config["use_nextalign"]:
    rule build_align:
        message:
            """
            Aligning sequences to {input.reference}
              - gaps relative to reference are considered real
            """
        input:
            sequences = rules.combine_samples.output.sequences,
            genemap = config["files"]["annotation"],
            reference = config["files"]["alignment_reference"]
        output:
            alignment = "builds/{build_name}/aligned.fasta",
            insertions = "builds/{build_name}/insertions.tsv",
            translations = expand("builds/{{build_name}}/translations/aligned.gene.{gene}.fasta", gene=config.get('genes', ['S']))
        params:
            outdir = "builds/{build_name}/translations",
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
            nextalign \
                --jobs={threads} \
                --reference {input.reference} \
                --genemap {input.genemap} \
                --genes {params.genes} \
                --sequences {input.sequences} \
                --output-dir {params.outdir} \
                --output-basename {params.basename} \
                --output-fasta {output.alignment} \
                --output-insertions {output.insertions} > {log} 2>&1
            """
else:
    rule build_align:
        message:
            """
            Aligning sequences from {input.sequences} to {input.reference}
            - gaps relative to reference are considered real
            """
        input:
            sequences = rules.combine_samples.output.sequences,
            reference = config["files"]["alignment_reference"]
        output:
            alignment = "builds/{build_name}/aligned.fasta"
        log:
            "logs/align_{build_name}.txt"
        benchmark:
            "benchmarks/align_{build_name}.txt"
        threads: 16
        conda: config["conda_environment"]
        shell:
            """
            mafft \
                --auto \
                --thread {threads} \
                --keeplength \
                --addfragments \
                {input.sequences} \
                {input.reference} > {output} 2> {log}
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
            lineages = "builds/{build_name}/pangolineages.csv",
        params:
            outdir = "builds/{build_name}",
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
            node_data = "builds/{build_name}/pangolineages.json"
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
        metadata = _get_unified_metadata
    output:
        metadata = "builds/{build_name}/metadata_adjusted.tsv"
    params:
        region = lambda wildcards: config["builds"][wildcards.build_name]["region"]
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
        tree = "builds/{build_name}/tree_raw.nwk"
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
        metadata = _get_metadata_by_wildcards
    output:
        tree = "builds/{build_name}/tree.nwk",
        node_data = "builds/{build_name}/branch_lengths.json"
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
        node_data = "builds/{build_name}/nt_muts.json"
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

rule haplotype_status:
    message: "Annotating haplotype status relative to {params.reference_node_name}"
    input:
        nt_muts = rules.ancestral.output.node_data
    output:
        node_data = "builds/{build_name}/haplotype_status.json"
    log:
        "logs/haplotype_status_{build_name}.txt"
    benchmark:
        "benchmarks/haplotype_status_{build_name}.txt"
    params:
        reference_node_name = config["reference_node_name"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/annotate-haplotype-status.py \
            --ancestral-sequences {input.nt_muts} \
            --reference-node-name {params.reference_node_name:q} \
            --output {output.node_data} 2>&1 | tee {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = config["files"]["reference"]
    output:
        node_data = "builds/{build_name}/aa_muts.json"
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
        node_data = "builds/{build_name}/aa_muts_explicit.json",
        translations = expand("builds/{{build_name}}/translations/aligned.gene.{gene}_withInternalNodes.fasta", gene=config.get('genes', ['S']))
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
if "use_nextalign" in config and config["use_nextalign"]:
    rule build_mutation_summary:
        message: "Summarizing {input.alignment}"
        input:
            alignment = rules.build_align.output.alignment,
            insertions = rules.build_align.output.insertions,
            translations = rules.build_align.output.translations,
            reference = config["files"]["alignment_reference"],
            genemap = config["files"]["annotation"]
        output:
            mutation_summary = "builds/{build_name}/mutation_summary.tsv"
        log:
            "logs/mutation_summary_{build_name}.txt"
        params:
            outdir = "builds/{build_name}/translations",
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
            alignments = "builds/{build_name}/translations/aligned.gene.S_withInternalNodes.fasta",
            distance_maps = ["defaults/distance_maps/S1.json"]
        params:
            genes = 'S',
            comparisons = ['root'],
            attribute_names = ['S1_mutations']
        output:
            node_data = "builds/{build_name}/distances.json"
        shell:
            """
            python scripts/mutation_counts.py \
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
        metadata = _get_metadata_by_wildcards
    output:
        node_data = "builds/{build_name}/traits.json"
    log:
        "logs/traits_{build_name}.txt"
    benchmark:
        "benchmarks/traits_{build_name}.txt"
    params:
        columns = _get_trait_columns_by_wildcards,
        sampling_bias_correction = _get_sampling_bias_correction_for_wildcards
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
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
        "builds/{build_name}/clades.tsv"
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
        clade_data = "builds/{build_name}/clades.json"
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
        clade_data = "builds/{build_name}/temp_emerging_lineages.json"
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
        clade_data = "builds/{build_name}/emerging_lineages.json"
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
        metadata = _get_metadata_by_wildcards
    output:
        colors = "builds/{build_name}/colors.tsv"
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
        metadata = _get_metadata_by_wildcards
    output:
        node_data = "builds/{build_name}/recency.json"
    log:
        "logs/recency_{build_name}.txt"
    benchmark:
        "benchmarks/recency_{build_name}.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
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
        metadata = _get_metadata_by_wildcards
    output:
        tip_frequencies_json = "builds/{build_name}/tip-frequencies.json"
    log:
        "logs/tip_frequencies_{build_name}.txt"
    benchmark:
        "benchmarks/tip_frequencies_{build_name}.txt"
    params:
        min_date = config["frequencies"]["min_date"],
        max_date = _get_max_date_for_frequencies,
        pivot_interval = config["frequencies"]["pivot_interval"],
        pivot_interval_units = config["frequencies"]["pivot_interval_units"],
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"]
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
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
        tree="builds/{build_name}/tree.nwk",
        frequencies="builds/{build_name}/tip-frequencies.json",
    output:
        node_data="builds/{build_name}/logistic_growth.json"
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

rule nucleotide_mutation_frequencies:
    message: "Estimate nucleotide mutation frequencies"
    input:
        alignment = rules.build_align.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = "builds/{build_name}/nucleotide_mutation_frequencies.json"
    log:
        "logs/nucleotide_mutation_frequencies_{build_name}.txt"
    benchmark:
        "benchmarks/nucleotide_mutation_frequencies_{build_name}.txt"
    params:
        min_date = config["frequencies"]["min_date"],
        minimal_frequency = config["frequencies"]["minimal_frequency"],
        pivot_interval = config["frequencies"]["pivot_interval"],
        stiffness = config["frequencies"]["stiffness"],
        inertia = config["frequencies"]["inertia"]
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        augur frequencies \
            --method diffusion \
            --alignments {input.alignment} \
            --gene-names nuc \
            --metadata {input.metadata} \
            --min-date {params.min_date} \
            --minimal-frequency {params.minimal_frequency} \
            --pivot-interval {params.pivot_interval} \
            --stiffness {params.stiffness} \
            --inertia {params.inertia} \
            --output {output.frequencies} 2>&1 | tee {log}
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
        rules.logistic_growth.output.node_data
    ]

    if "use_nextalign" in config and config["use_nextalign"]:
        inputs.append(rules.aa_muts_explicit.output.node_data)
        inputs.append(rules.distances.output.node_data)
    if "run_pangolin" in config and config["run_pangolin"]:
        inputs.append(rules.make_pangolin_node_data.output.node_data)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]

    return inputs

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_by_wildcards,
        auspice_config = lambda w: config["builds"][w.build_name]["auspice_config"] if "auspice_config" in config["builds"][w.build_name] else config["files"]["auspice_config"],
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"][w.build_name] else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"],
        description = lambda w: config["builds"][w.build_name]["description"] if "description" in config["builds"][w.build_name] else config["files"]["description"]
    output:
        auspice_json = "builds/{build_name}/ncov_with_accessions.json",
        root_sequence_json = "builds/{build_name}/ncov_with_accessions_root-sequence.json"
    log:
        "logs/export_{build_name}.txt"
    benchmark:
        "benchmarks/export_{build_name}.txt"
    params:
        title = export_title
    resources:
        # Memory use scales primarily with the size of the metadata file.
        mem_mb=lambda wildcards, input: 15 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        export AUGUR_RECURSION_LIMIT=10000 && augur export v2 \
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
        auspice_json = "builds/{build_name}/ncov_with_branch_labels.json"
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

rule incorporate_travel_history:
    message: "Adjusting main auspice JSON to take into account travel history"
    input:
        auspice_json = rules.add_branch_labels.output.auspice_json,
        colors = lambda w: config["builds"][w.build_name]["colors"] if "colors" in config["builds"][w.build_name] else ( config["files"]["colors"] if "colors" in config["files"] else rules.colors.output.colors.format(**w) ),
        lat_longs = config["files"]["lat_longs"]
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = "builds/{build_name}/ncov_with_accessions_and_travel_branches.json"
    log:
        "logs/incorporate_travel_history_{build_name}.txt"
    benchmark:
        "benchmarks/incorporate_travel_history_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule finalize:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = lambda w: rules.add_branch_labels.output.auspice_json if config.get("skip_travel_history_adjustment", False) else rules.incorporate_travel_history.output.auspice_json,
        frequencies = rules.tip_frequencies.output.tip_frequencies_json,
        root_sequence_json = rules.export.output.root_sequence_json
    output:
        auspice_json = "auspice/ncov_{build_name}.json",
        tip_frequency_json = "auspice/ncov_{build_name}_tip-frequencies.json",
        root_sequence_json = "auspice/ncov_{build_name}_root-sequence.json"
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
