# How to run: if no region is specified, it'll run a subsampled global build (120 per division)
# If a region is selected, it'll do 280/division for that region, and 20/division in the rest of the world
#       -- preferentially sequences near the focal sequences
#
# To run a regional build, be sure to update it to the list in Snakefile
#
# You can run all builds in parallel!
#   snakemake -s Snakefile_Regions --cores 2 all_regions
#
# Or you can specify final or intermediate output files like so:
#   snakemake -s Snakefile_Regions --cores 2 auspice/ncov_europe.json (subsampled regional focal)
#   snakemake -s Snakefile_Regions --cores 2 auspice/ncov.json (subsampled global)
#
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall -s Snakefile_Regions
#   snakemake -s Snakefile_Regions clean_export_regions
#   snakemake -s Snakefile_Regions export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake -s Snakefile_Regions
# to produce the final Auspice files!

rule all_regions:
    input:
        auspice_json = expand("auspice/ncov_{region}.json", region=REGIONS),
        tip_frequencies_json = expand("auspice/ncov_{region}_tip-frequencies.json", region=REGIONS),
        dated_auspice_json = expand("auspice/ncov_{region}_{date}.json", date=get_todays_date(), region=REGIONS),
        dated_tip_frequencies_json = expand("auspice/ncov_{region}_{date}_tip-frequencies.json", date=get_todays_date(), region=REGIONS),
        auspice_json_gisaid = expand("auspice/ncov_{region}_gisaid.json", region=REGIONS),
        auspice_json_zh = expand("auspice/ncov_{region}_zh.json", region=REGIONS)

# This cleans out files to allow re-run of 'normal' run (not ZH or GISAID)
# with `export` to check lat-longs & orderings
# However, *removes* the ZH & GISAID files so that when doing final run after all
# errors are cleared, these builds will also be rebuilt!
rule clean_export_regions:
    message: "Removing export files: {params}"
    params:
        "results/ncov_with_accessions*.json",
        "results/ncov_gisaid_with_accessions*.json",
        "results/ncov_zh_with_accessions*.json",
        "auspice/ncov*_gisaid.json",
        "auspice/ncov*_zh.json",
        "config/colors*.tsv"
    conda: "../envs/nextstrain.yaml"
    shell:
        "rm {params}"

# Allows 'normal' run of export to be forced to correct lat-long & ordering
# Just runs this, not ZH & GISAID, to speed up & reduce errors.
rule export_all_regions:
    input:
        colors_file = expand("config/colors_{region}.tsv", region=REGIONS),
        auspice_json = expand(REGION_PATH + "ncov_with_accessions.json", region=REGIONS),

rule subsample_focus:
    message:
        """
        Subsample all sequences into a focal set for {wildcards.region} with {params.sequences_per_group} per region
        """
    input:
        sequences = rules.mask.output.alignment,
        metadata = rules.download.output.metadata,
        include = config["files"]["include"]
    output:
        sequences = REGION_PATH + "subsample_focus.fasta"
    params:
        group_by = config["subsample_focus"]["group_by"],
        sequences_per_group = config["subsample_focus"]["seq_per_group_regional"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude-where region!={wildcards.region} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences} \
        """

rule make_priorities:
    message:
        """
        determine priority for inclusion in as phylogenetic context by
        genetic similiarity to sequences in focal set for region '{wildcards.region}'.
        """
    input:
        alignment = rules.mask.output.alignment,
        metadata = rules.download.output.metadata,
        focal_alignment = rules.subsample_focus.output.sequences
    output:
        priorities = REGION_PATH + "subsampling_priorities.tsv"
    resources:
        mem_mb = 4000
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/priorities.py --alignment {input.alignment} \
            --metadata {input.metadata} \
            --focal-alignment {input.focal_alignment} \
            --output {output.priorities}
        """

rule subsample_context:
    message:
        """
        Subsample the non-focal sequences to provide phylogenetic context for the region '{wildcards.region}' using {params.sequences_per_group} per {params.group_by}.
        """
    input:
        sequences = rules.mask.output.alignment,
        metadata = rules.download.output.metadata,
        priorities = rules.make_priorities.output.priorities
    output:
        sequences = REGION_PATH + "subsample_context.fasta"
    params:
        group_by = config["subsample_context"]["group_by"],
        sequences_per_group = config["subsample_context"]["sequences_per_group"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --exclude-where region={wildcards.region} \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --priority {input.priorities} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --output {output.sequences}
        """

rule subsample_regions:
    message:
        """
        Combine and deduplicate FASTAs
        """
    input:
        rules.subsample_focus.output.sequences,
        rules.subsample_context.output.sequences
    output:
        alignment = REGION_PATH + "subsampled_alignment.fasta"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/combine-and-dedup-fastas.py \
            --input {input} \
            --output {output}
        """

rule adjust_metadata_regions:
    message:
        """
        Adjusting metadata for region '{wildcards.region}'
        """
    input:
        metadata = rules.download.output.metadata
    output:
        metadata = REGION_PATH + "metadata_adjusted.tsv"
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        python3 scripts/adjust_regional_meta.py \
            --region "$region" \
            --metadata {input.metadata} \
            --output {output.metadata}
        """

#
# Deployment and error handlers, including Slack messaging integrations.
#

from os import environ

SLACK_TOKEN   = environ["SLACK_TOKEN"]   = config["slack_token"]   or ""
SLACK_CHANNEL = environ["SLACK_CHANNEL"] = config["slack_channel"] or ""

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
        *rules.all_regions.input
    params:
        slack_message = f"Deployed <https://nextstrain.org/staging/ncov|nextstrain.org/staging/ncov> {deploy_origin}",
        s3_staging_url = config["s3_staging_url"]
    conda: "../envs/nextstrain.yaml"
    shell:
        """
        nextstrain deploy {params.s3_staging_url:q} {input:q}

        if [[ -n "$SLACK_TOKEN" && -n "$SLACK_CHANNEL" ]]; then
            curl https://slack.com/api/chat.postMessage \
                --header "Authorization: Bearer $SLACK_TOKEN" \
                --form-string channel="$SLACK_CHANNEL" \
                --form-string text={params.slack_message:q} \
                --fail --silent --show-error \
                --include
        fi
        """

onerror:
    slack_message = f"Build {deploy_origin} failed."

    if SLACK_TOKEN and SLACK_CHANNEL:
        shell(f"""
            curl https://slack.com/api/files.upload \
                --header "Authorization: Bearer $SLACK_TOKEN" \
                --form-string channels="$SLACK_CHANNEL" \
                --form-string initial_comment={{slack_message:q}} \
                --form file=@{{log:q}} \
                --form filetype=text \
                --fail --silent --show-error \
                --include
        """)
