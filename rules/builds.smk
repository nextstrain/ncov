rule download:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = config["sequences"],
        curated_alignment = config["curated_alignment"],
        metadata = config["metadata"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv {output.metadata:q}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta {output.sequences:q}
        aws s3 cp s3://nextstrain-ncov-private/testing_curated_alignment.fasta {output.curated_alignment:q}
        """

rule filter:
    message:
        """
        TEMPORARY -- filtering to <100 genomes, the majority of which are going to be
        part of the hand-curated & maintained alignment. This is only to have a reproducible
        test set for development. 
        
        These sequences have been chosen from a monophyletic clade in the non-subsampled build which had
        a subclade which appeared to be poorly aligned. The small number allows me to view all the rows of
        the alignment on my screen
        
        The idea for the future is that this step would
        only filter out incompleted sequences -- i.e. date missing, too-many-Ns etc etc.
        The curated alignment would already have been filtered, and include as many sequences
        as possible (i.e. would include a large fraction of those output from this rule...)

        """
    input:
        sequences = rules.download.output.sequences,
        metadata = rules.download.output.metadata,
        include = config["files"]["include"],
    output:
        sequences = "results/filtered.fasta"
    params:
        min_length = 10000000
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --min-length {params.min_length} \
            --output {output.sequences}
        """


rule align:
    message:
        """
        Aligning filtered sequences to a hand curated alignment
          - gaps relative to reference are considered real
        """
    input:
        curated = rules.download.output.curated_alignment,
        sequences = rules.filter.output.sequences,
        reference = config["files"]["reference"]
    output:
        alignment = "results/all_aligned.fasta"
    benchmark:
        "benchmarks/align.txt"
    threads: 2
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur align \
            --existing-alignment {input.curated} \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference \
            --fill-gaps
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
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"],
        mask_sites = config["mask"]["mask_sites"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """

def _get_alignments_for_tree(wildcards):
    """Use all sequences for global builds. Use a specific subsampled set of
    sequences for regional builds.
    """
    if wildcards.region == "global":
        return rules.mask.output.alignment
    else:
        return rules.subsample_regions.output.alignment

rule tree:
    message: "Building tree"
    input:
        alignment = _get_alignments_for_tree
    output:
        tree = REGION_PATH + "tree_raw.nwk"
    benchmark:
        "benchmarks/tree_{region}.txt"
    threads: 4
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

def _get_metadata_by_wildcards(wildcards):
    if wildcards.region == "global":
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
        alignment = _get_alignments_for_tree,
        metadata = _get_metadata_by_wildcards
    output:
        tree = REGION_PATH + "tree.nwk",
        node_data = REGION_PATH + "branch_lengths.json"
    benchmark:
        "benchmarks/refine_{region}.txt"
    threads: 1
    params:
        root = config["refine"]["root"],
        clock_rate = config["refine"]["clock_rate"],
        clock_std_dev = config["refine"]["clock_std_dev"],
        coalescent = config["refine"]["coalescent"],
        date_inference = config["refine"]["date_inference"],
        divergence_unit = config["refine"]["divergence_unit"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"]
    conda: "../envs/nextstrain.yaml"
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
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - not inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output
    output:
        node_data = REGION_PATH + "nt_muts.json"
    params:
        inference = config["ancestral"]["inference"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous
        """

rule haplotype_status:
    message: "Annotating haplotype status relative to {params.reference_node_name}"
    input:
        nt_muts = rules.ancestral.output.node_data
    output:
        node_data = REGION_PATH + "haplotype_status.json"
    params:
        reference_node_name = config["reference_node_name"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/annotate-haplotype-status.py \
            --ancestral-sequences {input.nt_muts} \
            --reference-node-name {params.reference_node_name:q} \
            --output {output.node_data}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = config["files"]["reference"]
    output:
        node_data = REGION_PATH + "aa_muts.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

def _get_sampling_trait_for_wildcards(wildcards):
    mapping = {"north-america": "country", "oceania": "country"} # TODO: switch to "division"
    return mapping[wildcards.region] if wildcards.region in mapping else "country"

def _get_exposure_trait_for_wildcards(wildcards):
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
        node_data = REGION_PATH + "traits.json",
    params:
        columns = _get_exposure_trait_for_wildcards,
        sampling_bias_correction = config["traits"]["sampling_bias_correction"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} \
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["files"]["clades"]
    output:
        clade_data = REGION_PATH + "clades.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule colors:
    message: "Constructing colors file"
    input:
        ordering = config["files"]["ordering"],
        color_schemes = config["files"]["color_schemes"],
        metadata = _get_metadata_by_wildcards
    output:
        colors = "config/colors_{region}.tsv"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata}
        """

rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata = _get_metadata_by_wildcards
    output:
        REGION_PATH + "recency.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output}
        """

rule tip_frequencies:
    message: "Estimating censored KDE frequencies for tips"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        tip_frequencies_json = "auspice/ncov_{region}_tip-frequencies.json"
    params:
        min_date = config["frequencies"]["min_date"],
        pivot_interval = config["frequencies"]["pivot_interval"],
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"]
    conda: "../envs/nextstrain.yaml"
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
            --output {output.tip_frequencies_json}
        """

def export_title(wildcards):
    region = wildcards.region

    if not region:
        return "Genomic epidemiology of novel coronavirus"
    elif region == "global":
        return "Genomic epidemiology of novel coronavirus - Global subsampling"
    else:
        region_title = region.replace("-", " ").title()
        return f"Genomic epidemiology of novel coronavirus - {region_title}-focused subsampling"

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = config["files"]["auspice_config"],
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"],
        description = config["files"]["description"],
        clades = rules.clades.output.clade_data,
        recency = rules.recency.output
    output:
        auspice_json = REGION_PATH + "ncov_with_accessions.json"
    params:
        title = export_title
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} {input.clades} {input.recency} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --title {params.title:q} \
            --description {input.description} \
            --output {output.auspice_json}
        """

rule export_gisaid:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = config["files"]["auspice_config_gisaid"],
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"],
        description = config["files"]["description"],
        clades = rules.clades.output.clade_data,
        recency = rules.recency.output
    output:
        auspice_json = REGION_PATH + "ncov_gisaid_with_accessions.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} {input.clades} {input.recency} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --output {output.auspice_json}
        """

rule export_zh:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = config["files"]["auspice_config_zh"],
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"],
        description = config["files"]["description_zh"],
        clades = rules.clades.output.clade_data,
        recency = rules.recency.output
    output:
        auspice_json = REGION_PATH + "ncov_zh_with_accessions.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} {input.clades} {input.recency} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --output {output.auspice_json}
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
        auspice_json = REGION_PATH + "ncov_with_accessions_and_travel_branches.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json}
        """

rule incorporate_travel_history_gisaid:
    message: "Adjusting GISAID auspice JSON to take into account travel history"
    input:
        auspice_json = rules.export_gisaid.output.auspice_json,
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"]
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = REGION_PATH + "ncov_gisaid_with_accessions_and_travel_branches.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json}
        """

rule incorporate_travel_history_zh:
    message: "Adjusting ZH auspice JSON to take into account travel history"
    input:
        auspice_json = rules.export_zh.output.auspice_json,
        colors = rules.colors.output.colors,
        lat_longs = config["files"]["lat_longs"]
    params:
        sampling = _get_sampling_trait_for_wildcards,
        exposure = _get_exposure_trait_for_wildcards
    output:
        auspice_json = REGION_PATH + "ncov_zh_with_accessions_and_travel_branches.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json}
        """

rule fix_colorings:
    message: "Remove extraneous colorings for main build"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{region}.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """

rule fix_colorings_gisaid:
    message: "Remove extraneous colorings for the GISAID build"
    input:
        auspice_json = rules.incorporate_travel_history_gisaid.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{region}_gisaid.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """

rule fix_colorings_zh:
    message: "Remove extraneous colorings for the Chinese language build"
    input:
        auspice_json = rules.incorporate_travel_history_zh.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{region}_zh.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """

rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.fix_colorings.output.auspice_json,
        tip_frequencies_json = rules.tip_frequencies.output.tip_frequencies_json
    output:
        dated_auspice_json = "auspice/ncov_{region}_{date}.json",
        dated_tip_frequencies_json = "auspice/ncov_{region}_{date}_tip-frequencies.json"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        """
