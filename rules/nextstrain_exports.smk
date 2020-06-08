# How to run: if no region is specified, it'll run a subsampled global build (120 per division)
# If a region is selected, it'll do 280/division for that region, and 20/division in the rest of the world
#       -- preferentially sequences near the focal sequences
#
# To run a regional build, be sure to update the list of regions in `config/nextstrain_config.yaml`.
#
# You can run all builds in parallel!
#   snakemake --profile profiles/nextstrain all_regions
#
# Or you can specify final or intermediate output files like so:
#   snakemake --profile profiles/nextstrain auspice/ncov_europe.json (subsampled regional focal)
#   snakemake --profile profiles/nextstrain auspice/ncov_global.json (subsampled global)
#
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile profiles/nextstrain
#   snakemake --profile profiles/nextstrain clean_export_regions
#   snakemake --profile profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile profiles/nextstrain all_regions
# to produce the final Auspice files!

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

rule all_regions:
    input:
        auspice_json = expand("auspice/ncov_{build_name}.json", build_name=BUILD_NAMES),
        tip_frequencies_json = expand("auspice/ncov_{build_name}_tip-frequencies.json", build_name=BUILD_NAMES),
        dated_auspice_json = expand("auspice/ncov_{build_name}_{date}.json", build_name=BUILD_NAMES, date=get_todays_date()),
        dated_tip_frequencies_json = expand("auspice/ncov_{build_name}_{date}_tip-frequencies.json", build_name=BUILD_NAMES, date=get_todays_date()),
        auspice_json_gisaid = expand("auspice/ncov_{build_name}_gisaid.json", build_name=BUILD_NAMES),
        auspice_json_zh = expand("auspice/ncov_{build_name}_zh.json", build_name=BUILD_NAMES)

# This cleans out files to allow re-run of 'normal' run (not ZH or GISAID)
# with `export` to check lat-longs & orderings
# However, *removes* the ZH & GISAID files so that when doing final run after all
# errors are cleared, these builds will also be rebuilt!
rule clean_export_regions:
    message: "Removing export files: {input}"
    params:
        *expand("results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        *expand("results/{build_name}/ncov_gisaid_with_accessions.json", build_name=BUILD_NAMES),
        *expand("results/{build_name}/ncov_zh_with_accessions.json", build_name=BUILD_NAMES),
        *expand("results/{build_name}/colors.tsv", build_name=BUILD_NAMES),
        "auspice/ncov*_gisaid.json",
        "auspice/ncov*_zh.json"
    conda: config["conda_environment"]
    shell:
        "rm -f {params}"

# Allows 'normal' run of export to be forced to correct lat-long & ordering
# Just runs this, not ZH & GISAID, to speed up & reduce errors.
# Runs an additional script to give a list of locations that need colors and/or lat-longs
rule export_all_regions:
    input:
        auspice_json = expand("results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        lat_longs = config["files"]["lat_longs"],
        metadata = [_get_metadata_by_build_name(build_name).format(build_name=build_name)
                    for build_name in BUILD_NAMES],
        colors = expand("results/{build_name}/colors.tsv", build_name=BUILD_NAMES),
    conda: config["conda_environment"]
    shell:
        """
        {python:q} ./scripts/check_missing_locations.py \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --latlong {input.lat_longs}
        """

#
# Rules for custom auspice exports for the Nextstrain team.
#

rule export_gisaid:
    message: "Exporting GISAID-related data files for for auspice"
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
        auspice_json = "results/{build_name}/ncov_gisaid_with_accessions.json"
    log:
        "logs/export_gisaid_{build_name}.txt"
    conda: config["conda_environment"]
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
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule export_zh:
    message: "Exporting ZH-related data files for for auspice"
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
        auspice_json = "results/{build_name}/ncov_zh_with_accessions.json"
    log:
        "logs/export_zh_{build_name}.txt"
    conda: config["conda_environment"]
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
            --output {output.auspice_json} 2>&1 | tee {log}
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
        auspice_json = "results/{build_name}/ncov_gisaid_with_accessions_and_travel_branches.json"
    log:
        "logs/incorporate_travel_history_gisaid_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        {python:q} ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json} 2>&1 | tee {log}
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
        auspice_json = "results/{build_name}/ncov_zh_with_accessions_and_travel_branches.json"
    log:
        "logs/incorporate_travel_history_zh_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        {python:q} ./scripts/modify-tree-according-to-exposure.py \
            --input {input.auspice_json} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --sampling {params.sampling} \
            --exposure {params.exposure} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule fix_colorings_gisaid:
    message: "Remove extraneous colorings for the GISAID build"
    input:
        auspice_json = rules.incorporate_travel_history_gisaid.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{build_name}_gisaid.json"
    log:
        "logs/fix_colorings_gisaid_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        {python:q} scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule fix_colorings_zh:
    message: "Remove extraneous colorings for the Chinese language build"
    input:
        auspice_json = rules.incorporate_travel_history_zh.output.auspice_json
    output:
        auspice_json = "auspice/ncov_{build_name}_zh.json"
    log:
        "logs/fix_colorings_zh_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        {python:q} scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.finalize.output.auspice_json,
        tip_frequencies_json = rules.tip_frequencies.output.tip_frequencies_json
    output:
        dated_auspice_json = "auspice/ncov_{build_name}_{date}.json",
        dated_tip_frequencies_json = "auspice/ncov_{build_name}_{date}_tip-frequencies.json"
    conda: config["conda_environment"]
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
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
    conda: config["conda_environment"]
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
