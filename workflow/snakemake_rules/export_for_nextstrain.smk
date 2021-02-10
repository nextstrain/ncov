# How to run: if no region is specified, it'll run a subsampled global build (120 per division)
# If a region is selected, it'll do 280/division for that region, and 20/division in the rest of the world
#       -- preferentially sequences near the focal sequences
#
# To run a regional build, be sure to update the list of regions in `config/nextstrain_profiles.yaml`.
#
# You can run all builds in parallel!
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
#
# Or you can specify final or intermediate output files like so:
#   snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_europe.json (subsampled regional focal)
#   snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_global.json (subsampled global)
#
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
#   snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
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
        dated_tip_frequencies_json = expand("auspice/ncov_{build_name}_{date}_tip-frequencies.json", build_name=BUILD_NAMES, date=get_todays_date())

# This cleans out files to allow re-run of 'normal' run with `export`
# to check lat-longs & orderings
rule clean_export_regions:
    message: "Removing export files: {input}"
    params:
        *expand("results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        *expand("results/{build_name}/colors.tsv", build_name=BUILD_NAMES)
    conda: config["conda_environment"]
    shell:
        "rm -f {params}"

# Allows 'normal' run of export to be forced to correct lat-long & ordering
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
        python3 ./scripts/check_missing_locations.py \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --latlong {input.lat_longs}
        """

rule all_mutation_frequencies:
    input: expand("results/{build_name}/nucleotide_mutation_frequencies.json", build_name=BUILD_NAMES)

#
# Rules for custom auspice exports for the Nextstrain team.
#

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


rule upload:
    message: "Uploading intermediate files for specified origins to {params.s3_bucket}"
    input:
        expand("results/aligned_{origin}.fasta", origin=config["S3_DST_ORIGINS"]),              # from `rule align`
        expand("results/sequence-diagnostics_{origin}.tsv", origin=config["S3_DST_ORIGINS"]),   # from `rule diagnostic`
        expand("results/flagged-sequences_{origin}.tsv", origin=config["S3_DST_ORIGINS"]),      # from `rule diagnostic`
        expand("results/to-exclude_{origin}.txt", origin=config["S3_DST_ORIGINS"]),             # from `rule diagnostic`
        expand("results/aligned-filtered_{origin}.fasta", origin=config["S3_DST_ORIGINS"]),     # from `rule refilter`
        expand("results/masked_{origin}.fasta", origin=config["S3_DST_ORIGINS"]),               # from `rule mask`
        expand("results/filtered_{origin}.fasta", origin=config["S3_DST_ORIGINS"]),             # from `rule filter`
    params:
        s3_bucket = config["S3_DST_BUCKET"],
        compression = config["S3_DST_COMPRESSION"]
    log:
        "logs/upload_gisaid.txt"
    run:
        for fname in input:
            cmd = f"./scripts/upload-to-s3 {fname} s3://{params.s3_bucket}/{os.path.basename(fname)}.{params.compression} | tee -a {log}"
            print("upload command:", cmd)
            shell(cmd)

onstart:
    slack_message = f"Build {deploy_origin} started."

    if SLACK_TOKEN and SLACK_CHANNEL:
        shell(f"""
            curl https://slack.com/api/chat.postMessage \
                --header "Authorization: Bearer $SLACK_TOKEN" \
                --form-string channel="$SLACK_CHANNEL" \
                --form-string text={{slack_message:q}} \
                --fail --silent --show-error \
                --include
        """)

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
