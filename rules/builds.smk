rule download:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = config["sequences"],
        metadata = config["metadata"]
    conda: config["conda_environment"]
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq >{output.metadata:q}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta.gz - | gunzip -cq > {output.sequences:q}
        """

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length}
        """
    input:
        sequences = rules.download.output.sequences,
        metadata = rules.download.output.metadata,
        include = config["files"]["include"],
        exclude = config["files"]["exclude"]
    output:
        sequences = "results/filtered.fasta"
    log:
        "logs/filtered.txt"
    params:
        min_length = config["filter"]["min_length"],
        exclude_where = config["filter"]["exclude_where"],
        group_by = config["filter"]["group_by"],
        sequences_per_group = config["filter"]["sequences_per_group"]
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where}\
            --min-length {params.min_length} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences} 2>&1 | tee {log}
        """

checkpoint partition_sequences:
    input:
        sequences = rules.filter.output.sequences
    output:
        split_sequences = directory("results/split_sequences/")
    log:
        "logs/partition_sequences.txt"
    params:
        sequences_per_group = config["partition_sequences"]["sequences_per_group"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences} 2>&1 | tee {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        Cluster:  {wildcards.cluster}
        """
    input:
        sequences = "results/split_sequences/{cluster}.fasta",
        reference = config["files"]["reference"]
    output:
        alignment = "results/split_alignments/{cluster}.fasta"
    log:
        "logs/align_{cluster}.txt"
    benchmark:
        "benchmarks/align_{cluster}.txt"
    threads: 2
    conda: config["conda_environment"]
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference \
            --fill-gaps 2>&1 | tee {log}
        """

def _get_alignments(wildcards):
    checkpoint_output = checkpoints.partition_sequences.get(**wildcards).output[0]
    return expand("results/split_alignments/{i}.fasta",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule aggregate_alignments:
    message: "Collecting alignments"
    input:
        alignments = _get_alignments
    output:
        alignment = "results/aligned.fasta"
    log:
        "logs/aggregate_alignments.txt"
    conda: config["conda_environment"]
    shell:
        """
        cat {input.alignments} > {output.alignment} 2> {log}
        """

rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.aggregate_alignments.output.alignment
    output:
        alignment = "results/masked.fasta"
    log:
        "logs/mask.txt"
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
            --output {output.alignment} 2>&1 | tee {log}
        """

def _get_subsampling_settings(wildcards):
    # Allow users to override default subsampling with their own settings keyed
    # by location type and name. For example, "region_europe" or
    # "country_iceland". Otherwise, default to settings for the location type.
    custom_subsampling_key = f"{wildcards.location_type}_{wildcards.location_name}"

    if custom_subsampling_key in config["subsampling"]:
        subsampling_settings = config["subsampling"][custom_subsampling_key]
    else:
        subsampling_settings = config["subsampling"][wildcards.location_type]

    if hasattr(wildcards, "subsample"):
        return subsampling_settings[wildcards.subsample]
    else:
        return subsampling_settings

def get_priorities(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return f"results/{wildcards.location_type}/{wildcards.location_name}/proximity_{subsampling_settings['priorities']['focus']}.tsv"
    else:
        # TODO: find a way to make the list of input files depend on config
        return config["files"]["include"]

def get_priority_argument(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    if "priorities" in subsampling_settings and subsampling_settings["priorities"]["type"] == "proximity":
        return "--priority " + get_priorities(wildcards)
    else:
        return ""

def _get_subsample_group_by(wildcards):
    return _get_subsampling_settings(wildcards)["group_by"]

def _get_specific_subsampling_setting(setting, optional=False):
    def _get_setting(wildcards):
        if optional:
            value = _get_subsampling_settings(wildcards).get(setting, "")
        else:
            value = _get_subsampling_settings(wildcards)[setting]

        if isinstance(value, str):
            value = value.format(**wildcards)
        else:
            return value

        if "geo_hierarchy" in config\
             and wildcards.location_type in config["geo_hierarchy"] \
             and wildcards.location_name in config["geo_hierarchy"][wildcards.location_type]:
            return value.format(**config["geo_hierarchy"][wildcards.location_type][wildcards.location_name])
        else:
            return value

    return _get_setting

rule subsample:
    message:
        """
        Subsample all sequences into a {wildcards.subsample} set for {wildcards.location_type} {wildcards.location_name} with {params.sequences_per_group} per {params.group_by}
        """
    input:
        sequences = rules.mask.output.alignment,
        metadata = rules.download.output.metadata,
        include = config["files"]["include"],
        priorities = get_priorities
    output:
        sequences = "results/{location_type}/{location_name}/sample-{subsample}.fasta"
    params:
        group_by = _get_specific_subsampling_setting("group_by"),
        sequences_per_group = _get_specific_subsampling_setting("seq_per_group"),
        exclude_argument = _get_specific_subsampling_setting("exclude", optional=True),
        include_argument = _get_specific_subsampling_setting("include", optional=True),
        priority_argument = get_priority_argument
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            {params.exclude_argument} \
            {params.include_argument} \
            {params.priority_argument} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule proximity_score:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for {wildcards.location_type} '{wildcards.location_name}'.
        """
    input:
        alignment = rules.mask.output.alignment,
        metadata = rules.download.output.metadata,
        focal_alignment = "results/{location_type}/{location_name}/sample-{focus}.fasta"
    output:
        priorities = "results/{location_type}/{location_name}/proximity_{focus}.tsv"
    log:
        "logs/subsampling_priorities_{location_type}_{location_name}_{focus}.txt"
    resources:
        mem_mb = 4000
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/priorities.py --alignment {input.alignment} \
            --metadata {input.metadata} \
            --focal-alignment {input.focal_alignment} \
            --output {output.priorities} 2>&1 | tee {log}
        """

def _get_subsampled_files(wildcards):
    subsampling_settings = _get_subsampling_settings(wildcards)

    return [
        f"results/{wildcards.location_type}/{wildcards.location_name}/sample-{subsample}.fasta"
        for subsample in subsampling_settings
    ]

rule combine_samples:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        _get_subsampled_files
    output:
        alignment = "results/{location_type}/{location_name}/subsampled_alignment.fasta"
    log:
        "logs/subsample_regions_{location_type}_{location_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py \
            --input {input} \
            --output {output} 2>&1 | tee {log}
        """

rule adjust_metadata_regions:
    message:
        """
        Adjusting metadata for {wildcards.location_type} '{wildcards.location_name}'
        """
    input:
        metadata = rules.download.output.metadata
    output:
        metadata = "results/{location_type}/{location_name}/metadata_adjusted.tsv"
    log:
        "logs/adjust_metadata_regions_{location_type}_{location_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/adjust_regional_meta.py \
            --{wildcards.location_type} "{wildcards.location_name}" \
            --metadata {input.metadata} \
            --output {output.metadata} 2>&1 | tee {log}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.combine_samples.output.alignment
    output:
        tree = "results/{location_type}/{location_name}/tree_raw.nwk"
    log:
        "logs/tree_{location_type}_{location_name}.txt"
    benchmark:
        "benchmarks/tree_{location_type}_{location_name}.txt"
    threads: 16
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
            --output {output.tree} \
            --nthreads {threads} 2>&1 | tee {log}
        """

def _get_metadata_by_wildcards(wildcards):
    if wildcards.location_name == "global":
        return rules.download.output.metadata
    else:
        return rules.adjust_metadata_regions.output.metadata

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
        alignment = rules.combine_samples.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        tree = "results/{location_type}/{location_name}/tree.nwk",
        node_data = "results/{location_type}/{location_name}/branch_lengths.json"
    log:
        "logs/refine_{location_type}_{location_name}.txt"
    benchmark:
        "benchmarks/refine_{location_type}_{location_name}.txt"
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
        clock_filter_iqd = config["refine"]["clock_filter_iqd"]
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
            --timetree \
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
        alignment = rules.combine_samples.output.alignment
    output:
        node_data = "results/{location_type}/{location_name}/nt_muts.json"
    log:
        "logs/ancestral_{location_type}_{location_name}.txt"
    params:
        inference = config["ancestral"]["inference"]
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
        node_data = "results/{location_type}/{location_name}/haplotype_status.json"
    log:
        "logs/haplotype_status_{location_type}_{location_name}.txt"
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
        node_data = "results/{location_type}/{location_name}/aa_muts.json"
    log:
        "logs/translate_{location_type}_{location_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

def _get_sampling_trait_for_wildcards(wildcards):
    # TODO: fix this for locations
    return "country"
    mapping = {"north-america": "country", "oceania": "country"} # TODO: switch to "division"
    return mapping[wildcards.region] if wildcards.region in mapping else "country"

def _get_exposure_trait_for_wildcards(wildcards):
    # TODO: fix this for locations
    return "country_exposure"
    mapping = {"north-america": "country_exposure", "oceania": "country_exposure"} # TODO: switch to "division_exposure"
    return mapping[wildcards.region] if wildcards.region in mapping else "country_exposure"

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
        node_data = "results/{location_type}/{location_name}/traits.json"
    log:
        "logs/traits_{location_type}_{location_name}.txt"
    params:
        columns = lambda w: config["traits"][w.region]["columns"],
        sampling_bias_correction = lambda w: config["traits"][w.region]["sampling_bias_correction"]
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

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["files"]["clades"]
    output:
        clade_data = "results/{location_type}/{location_name}/clades.json"
    log:
        "logs/clades_{location_type}_{location_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data} 2>&1 | tee {log}
        """

rule colors:
    message: "Constructing colors file"
    input:
        ordering = config["files"]["ordering"],
        color_schemes = config["files"]["color_schemes"],
        metadata = _get_metadata_by_wildcards
    output:
        colors = "config/colors_{location_type}_{location_name}.tsv"
    log:
        "logs/colors_{location_type}_{location_name}.txt"
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
        node_data = "results/{location_type}/{location_name}/recency.json"
    log:
        "logs/recency_{location_type}_{location_name}.txt"
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
        tip_frequencies_json = "auspice/ncov_{location_type}_{location_name}_tip-frequencies.json"
    log:
        "logs/tip_frequencies_{location_type}_{location_name}.txt"
    params:
        min_date = config["frequencies"]["min_date"],
        pivot_interval = config["frequencies"]["pivot_interval"],
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"]
    conda: config["conda_environment"]
    shell:
        """
        augur frequencies \
            --method kde \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --min-date {params.min_date} \
            --pivot-interval {params.pivot_interval} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --output {output.tip_frequencies_json} 2>&1 | tee {log}
        """

def export_title(wildcards):
    location_name = wildcards.location_name

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
        rules.clades.output.clade_data,
        rules.recency.output.node_data
    ]

    # TODO: reenable traits
    #if wildcards.region in config["traits"]:
    #    inputs.append(rules.traits.output.node_data)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_by_wildcards,
        auspice_config = config["files"]["auspice_config"],
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"],
        description = config["files"]["description"]
    output:
        auspice_json = "results/{location_type}/{location_name}/ncov_with_accessions.json"
    log:
        "logs/export_{location_type}_{location_name}.txt"
    params:
        title = export_title
    conda: config["conda_environment"]
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule incorporate_travel_history:
    message: "Adjusting main auspice JSON to take into account travel history"
    input:
        auspice_json = rules.export.output.auspice_json,
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"]
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = "results/{location_type}/{location_name}/ncov_with_accessions_and_travel_branches.json"
    log:
        "logs/incorporate_travel_history_{location_type}_{location_name}.txt"
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

rule fix_colorings:
    message: "Remove extraneous colorings for main build"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{location_type}_{location_name}.json"
    log:
        "logs/fix_colorings_{location_type}_{location_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """
