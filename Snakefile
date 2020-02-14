import os
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider()
BUCKET = "nextstrain-ncov-private/"

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

rule all:
    input:
        auspice_json = "auspice/ncov.json",
        dated_auspice_json = expand("auspice/ncov_{date}.json", date=get_todays_date()),
        auspice_json_gisaid = "auspice/ncov_gisaid.json",
        auspice_json_zh = "auspice/ncov_zh.json"

rule files:
    params:
        input_fasta = "data/ncov.fasta",
        include = "config/include.txt",
        exclude = "config/exclude.txt",
        reference = "config/reference.gb",
        outgroup = "config/outgroup.fasta",
        weights = "config/weights.tsv",
        auspice_config = "config/auspice_config.json",
        auspice_config_gisaid = "config/auspice_config_gisaid.json",
        auspice_config_zh = "config/auspice_config_zh.json",
        colors = "config/colors.tsv",
        lat_longs = "config/lat_longs.tsv",
        description = "config/description.md",
        description_zh = "config/description_zh.md"

files = rules.files.params

rule download:
    message: "Downloading metadata and fasta files from S3"
    input:
        metadata = S3.remote(BUCKET + "metadata.tsv"),
        sequences = S3.remote(BUCKET + "sequences.fasta")
    output:
        metadata = "data/metadata.tsv",
        sequences = "data/sequences.fasta"
    shell:
        """
        cp {input.metadata:q} {output.metadata:q}
        cp {input.sequences:q} {output.sequences:q}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
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
        group_by = "country",
        sequences_per_group = 100,
        min_length = 5000,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
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
        mask_from_beginning = 130,
        mask_from_end = 15,
        mask_sites = 18529
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
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
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
        root = "Wuhan/WIV04/2019 Wuhan/WIV06/2019",
        clock_rate = 0.0005,
        clock_std_dev = 0.0003,
        coalescent = "skyline",
        date_inference = "marginal",
        divergence_unit = "mutations"
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
            --no-covariance
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
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
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
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
        tree = rules.refine.output.tree,
        metadata = rules.download.output.metadata,
        weights = files.weights
    output:
        node_data = "results/traits.json",
    params:
        columns = "division",
        sampling_bias_correction = 2
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
        colors = files.colors,
        lat_longs = files.lat_longs,
        description = files.description
    output:
        auspice_json = "results/ncov_with_accessions.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} \
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
        auspice_config = files.auspice_config_gisaid,
        colors = files.colors,
        lat_longs = files.lat_longs,
        description = files.description
    output:
        auspice_json = "results/ncov_gisaid_with_accessions.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
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
        auspice_config = files.auspice_config_zh,
        colors = files.colors,
        lat_longs = files.lat_longs,
        description = files.description_zh
    output:
        auspice_json = "results/ncov_zh_with_accessions.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --output {output.auspice_json}
        """

rule fix_colorings:
    message: "Remove extraneous colorings"
    input:
        auspice_json = rules.export.output.auspice_json
    output:
        auspice_json = "auspice/ncov.json"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """

rule fix_colorings_gisaid:
    message: "Remove extraneous colorings"
    input:
        auspice_json = rules.export_gisaid.output.auspice_json
    output:
        auspice_json = "auspice/ncov_gisaid.json"
    shell:
        """
        python scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json}
        """

rule fix_colorings_zh:
    message: "Remove extraneous colorings"
    input:
        auspice_json = rules.export_zh.output.auspice_json
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
        auspice_json = rules.export.output.auspice_json
    output:
        dated_auspice_json = rules.all.input.dated_auspice_json
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
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


rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
