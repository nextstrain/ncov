from os import environ
from socket import getfqdn
from getpass import getuser

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

configfile: "config/Snakefile.yaml"

rule all:
    input:
        auspice_json = "auspice/ncov.json",
        tip_frequencies_json = "auspice/ncov_tip-frequencies.json",
        dated_auspice_json = expand("auspice/ncov_{date}.json", date=get_todays_date()),
        dated_tip_frequencies_json = expand("auspice/ncov_{date}_tip-frequencies.json", date=get_todays_date()),
        auspice_json_gisaid = "auspice/ncov_gisaid.json",
        auspice_json_zh = "auspice/ncov_zh.json"

rule files:
    params:
        include = "config/include.txt",
        exclude = "config/exclude.txt",
        reference = "config/reference.gb",
        outgroup = "config/outgroup.fasta",
        weights = "config/weights.tsv",
        ordering = "config/ordering.tsv",
        color_schemes = "config/color_schemes.tsv",
        auspice_config = "config/auspice_config.json",
        auspice_config_gisaid = "config/auspice_config_gisaid.json",
        auspice_config_zh = "config/auspice_config_zh.json",
        lat_longs = "config/lat_longs.tsv",
        description = "config/description.md",
        description_zh = "config/description_zh.md",
        clades = "config/clades.tsv"

files = rules.files.params

rule download:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = config["sequences"],
        metadata = config["metadata"]
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta {output.sequences:q}
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv {output.metadata:q}
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
        include = files.include,
        exclude = files.exclude
    output:
        sequences = "results/filtered.fasta"
    params:
        min_length = 25000,
        group_by = "division year month",
        sequences_per_group = 300,
        exclude_where = "date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04'"
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
            --output {output.sequences}
        """

checkpoint partition_sequences:
    input:
        sequences = rules.filter.output.sequences
    output:
        split_sequences = directory("results/split_sequences/pre/")
    params:
        sequences_per_group = 150
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences}
        """

rule partitions_intermediate:
    message:
        """
        partitions_intermediate: Copying sequence fastas
        {wildcards.cluster}
        """
    input:
        "results/split_sequences/pre/{cluster}.fasta"
    output:
        "results/split_sequences/post/{cluster}.fasta"
    shell:
        "cp {input} {output}"

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        {wildcards.cluster}
        """
    input:
        sequences = rules.partitions_intermediate.output,
        reference = files.reference
    output:
        alignment = "results/split_alignments/{cluster}.fasta"
    threads: 2
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference \
            --fill-gaps
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
    shell:
        """
        cat {input.alignments} > {output.alignment}
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
    params:
        mask_from_beginning = 130,
        mask_from_end = 50,
        mask_sites = "18529 29849 29851 29853"
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
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
        alignment = rules.mask.output,
        metadata = rules.download.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root = "Wuhan-Hu-1/2019 Wuhan/WH01/2019",
        clock_rate = 0.0008,
        clock_std_dev = 0.0004,
        coalescent = "skyline",
        date_inference = "marginal",
        divergence_unit = "mutations",
        clock_filter_iqd = 4
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
        tree = "results/tree.nwk",
        alignment = rules.mask.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
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
        node_data = "results/haplotype_status.json"
    params:
        reference_node_name = "USA/WA1/2020"
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
        tree = "results/tree.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = "results/tree.nwk",
        metadata = rules.download.output.metadata,
        weights = files.weights
    output:
        node_data = "results/traits.json",
    params:
        columns = "division_exposure",
        sampling_bias_correction = 2.5
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --weights {input.weights} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction} \
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = "results/tree.nwk",
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        clade_data = "results/clades.json"
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
        ordering = files.ordering,
        color_schemes = files.color_schemes
    output:
        colors = "config/colors.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}
        """

rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata = rules.download.output.metadata
    output:
        "results/recency.json"
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
        metadata = rules.download.output.metadata
    output:
        tip_frequencies_json = "auspice/ncov_tip-frequencies.json"
    params:
        min_date = 2020.0,
        pivot_interval = 1,
        narrow_bandwidth = 0.05,
        proportion_wide = 0.0
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

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.download.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = files.auspice_config,
        colors = rules.colors.output.colors,
        lat_longs = files.lat_longs,
        description = files.description,
        clades = "results/clades.json",
        recency = rules.recency.output
    output:
        auspice_json = "results/ncov_with_accessions.json"
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

rule export_gisaid:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.download.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = files.auspice_config_gisaid,
        colors = rules.colors.output.colors,
        lat_longs = files.lat_longs,
        description = files.description,
        clades = rules.clades.output.clade_data,
        recency = rules.recency.output
    output:
        auspice_json = "results/ncov_gisaid_with_accessions.json"
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
        metadata = rules.download.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        traits = rules.traits.output.node_data,
        auspice_config = files.auspice_config_zh,
        colors = rules.colors.output.colors,
        lat_longs = files.lat_longs,
        description = files.description_zh,
        clades = rules.clades.output.clade_data,
        recency = rules.recency.output
    output:
        auspice_json = "results/ncov_zh_with_accessions.json"
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
        lat_longs = files.lat_longs,
        colors = rules.colors.output.colors,
        auspice_json = rules.export.output.auspice_json
    output:
        auspice_json = "results/ncov_with_accessions_and_travel_branches.json"
    shell:
        """
        python3 ./scripts/modify_tree_according_to_division_exposure.py \
            {input.auspice_json} {input.colors} {input.lat_longs} {output.auspice_json}
        """

rule incorporate_travel_history_gisaid:
    message: "Adjusting GISAID auspice JSON to take into account travel history"
    input:
        lat_longs = files.lat_longs,
        colors = rules.colors.output.colors,
        auspice_json = rules.export_gisaid.output.auspice_json
    output:
        auspice_json = "results/ncov_gisaid_with_accessions_and_travel_branches.json"
    shell:
        """
        python3 ./scripts/modify_tree_according_to_division_exposure.py \
            {input.auspice_json} {input.colors} {input.lat_longs} {output.auspice_json}
        """

rule incorporate_travel_history_zh:
    message: "Adjusting ZH auspice JSON to take into account travel history"
    input:
        lat_longs = files.lat_longs,
        colors = rules.colors.output.colors,
        auspice_json = rules.export_zh.output.auspice_json
    output:
        auspice_json = "results/ncov_zh_with_accessions_and_travel_branches.json"
    shell:
        """
        python3 ./scripts/modify_tree_according_to_division_exposure.py \
            {input.auspice_json} {input.colors} {input.lat_longs} {output.auspice_json}
        """

rule fix_colorings:
    message: "Remove extraneous colorings for main build"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json
    output:
        auspice_json = "auspice/ncov.json"
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
        auspice_json = "auspice/ncov_gisaid.json"
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
        auspice_json = "auspice/ncov_zh.json"
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
        dated_auspice_json = rules.all.input.dated_auspice_json,
        dated_tip_frequencies_json = rules.all.input.dated_tip_frequencies_json
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        """

rule poisson_tmrca:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.download.output.metadata,
        nt_muts = rules.ancestral.output.node_data
    output:
        "figures/ncov_poisson-tmrca.png"
    shell:
        """
        python scripts/tmrca_estimate.py --tree {input.tree} --metadata {input.metadata} --node-data {input.nt_muts} --output {output}
        """

rule branching_process_R0:
    params:
        infectious_period = 10, # days
        population = [6000, 30000, 150000],
        start_recent = "2019-12-01",
        start_early = "2019-11-01"
    output:
        "figures/ncov_branching-R0-recent.png",
        "figures/ncov_branching-R0-early.png"
    shell:
        """
        python scripts/branching_process.py --infectious-period {params.infectious_period}\
                    --start {params.start_recent} \
                    --population {params.population} \
                    --output {output[0]} &&\
        python scripts/branching_process.py --infectious-period {params.infectious_period}\
                    --start {params.start_early} \
                    --population {params.population} \
                    --output {output[1]}
        """

try:
    deploy_origin = (
        f"from AWS Batch job `{environ['AWS_BATCH_JOB_ID']}`"
        if environ.get("AWS_BATCH_JOB_ID") else
        f"by the hands of {getuser()}@{getfqdn()}"
    )
except:
    # getuser() and getfqdn() may not always succeed, and this catch-all except
    # means that the Snakefile won't crash.
    deploy_origin = "by an unknown identity"

rule deploy_to_staging:
    input:
        *rules.all.input
    params:
        slack_message = json.dumps({"text":f"Deployed <https://nextstrain.org/staging/ncov|nextstrain.org/staging/ncov> {deploy_origin}"}),
        slack_webhook = config["slack_webhook"] or "",
        s3_staging_url = config["s3_staging_url"]
    shell:
        """
        nextstrain deploy {params.s3_staging_url:q} {input:q}

        if [[ -n "{params.slack_webhook}" ]]; then
            curl {params.slack_webhook:q} \
                --header 'Content-type: application/json' \
                --data-raw {params.slack_message:q} \
                --fail --silent --show-error \
                --include
        fi
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
