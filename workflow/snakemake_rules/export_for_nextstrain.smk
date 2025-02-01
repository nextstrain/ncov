# How to run: if no region is specified, it'll run a subsampled global build (120 per division)
# If a region is selected, it'll do 280/division for that region, and 20/division in the rest of the world
#       -- preferentially sequences near the focal sequences
#
# To run a regional build, be sure to update the list of regions in `config/nextstrain_profiles.yaml`.
#
# You can run all builds in parallel!
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid all_regions
#
# Or you can specify final or intermediate output files like so:
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid auspice/ncov_europe.json (subsampled regional focal)
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid auspice/ncov_global.json (subsampled global)
#
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain-gisaid
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain-gisaid all_regions
# to produce the final Auspice files!

import re
import requests
import json
from math import ceil
from workflow.lib.persistent_dict import PersistentDict, NoSuchEntryError

ruleorder: dated_json > finalize

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

rule all_regions:
    input:
        auspice_json = expand("auspice/{prefix}_{build_name}.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES),
        tip_frequencies_json = expand("auspice/{prefix}_{build_name}_tip-frequencies.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES),
        root_sequence_json = expand("auspice/{prefix}_{build_name}_root-sequence.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES),
        dated_auspice_json = expand("auspice/{prefix}_{build_name}_{date}.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES, date=config.get("build_date", get_todays_date())),
        dated_tip_frequencies_json = expand("auspice/{prefix}_{build_name}_{date}_tip-frequencies.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES, date=config.get("build_date", get_todays_date())),
        dated_root_sequence_json = expand("auspice/{prefix}_{build_name}_{date}_root-sequence.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES, date=config.get("build_date", get_todays_date()))

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

# Build specific metadata
rule extract_meta:
    input:
        alignment = rules.build_align.output.alignment,
        metadata="results/{build_name}/metadata_adjusted.tsv.xz",
    output:
        metadata = "results/{build_name}/extracted_metadata.tsv"
    run:
        from Bio import SeqIO
        import pandas as pd

        seq_names = [s.id for s in SeqIO.parse(input.alignment, 'fasta')]
        all_meta = pd.read_csv(input.metadata, sep='\t', index_col=0, dtype=str)
        extracted_meta = all_meta.loc[seq_names]
        extracted_meta.to_csv(output.metadata, sep='\t')


# Allows 'normal' run of export to be forced to correct lat-long & ordering
# Runs an additional script to give a list of locations that need colors and/or lat-longs
rule export_all_regions:
    input:
        auspice_json = expand("results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        lat_longs = config["files"]["lat_longs"],
        metadata=expand("results/{build_name}/metadata_adjusted.tsv.xz", build_name=BUILD_NAMES),
        colors = expand("results/{build_name}/colors.tsv", build_name=BUILD_NAMES),
    benchmark:
        "benchmarks/export_all_regions.txt"
    resources:
        # Memory use scales primarily with the size of the metadata file.
        # Compared to other rules, this rule loads metadata as a pandas
        # DataFrame instead of a dictionary, so it uses much less memory.
        mem_mb=lambda wildcards, input: ceil(5 * (input.metadata.size / 1024 / 1024))
    conda: config["conda_environment"]
    shell:
        r"""
        python3 ./scripts/check_missing_locations.py \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --latlong {input.lat_longs}
        """


#
# Rule for generating a per-build auspice config
#
rule auspice_config:
    """
    This rule is only intended to be run with `nextstrain-open` or `nextstrain-gisaid`
    profiles!
    """
    message: "Making a custom auspice config."
    output:
        "results/{build_name}/auspice_config.json"
    benchmark:
        "benchmarks/make_auspice_config_{build_name}.txt"
    run:
        input_set = set(config['inputs']) - {"references"}
        build_name = wildcards.build_name

        if "_" in build_name:
            build_region, build_timespan = build_name.split("_")
        else:
            build_region = build_name
            build_timespan = ""

        ## What are the parameters which vary across builds?
        ## Note: set a value to `None` to exclude it from the produced config JSON
        default_geo_resolution = "division" if build_region in ["north-america", "oceania"] else "country"
        default_map_triplicate = True if build_region in ["reference", "global"] else False
        if input_set == {"gisaid"}:
            data_provenance = [{"name": "GISAID"}]
            gisaid_epi_isl_coloring = {"key": "gisaid_epi_isl", "type": "categorical"}
            location_coloring = {"key": "location", "title": "Location", "type": "categorical"}
            location_filter = "location"
            originating_lab_filter = "originating_lab"
            submitting_lab_filter  = "submitting_lab"
        elif input_set == {"open"}:
            data_provenance = [{"name": "GenBank", "url": "https://www.ncbi.nlm.nih.gov/genbank/"}]
            gisaid_epi_isl_coloring = None
            location_coloring = None
            location_filter = None
            originating_lab_filter = None
            submitting_lab_filter  = None
        else:
            raise Exception(f"rule auspice_config doesn't know how to handle inputs: {input_set}")

        data = {
            "build_url": "https://github.com/nextstrain/ncov",
            "maintainers": [
                {
                "name": "the Nextstrain team",
                "url": "https://nextstrain.org/"
                }
            ],
            "data_provenance": data_provenance,
            "colorings": [
                {
                    "key": "emerging_lineage",
                    "title": "Emerging Lineage",
                    "type": "categorical"
                },
                {
                    "key": "pango_lineage",
                    "title": "Pangolin Pango Lineage",
                    "type": "categorical"
                },
                {
                    "key": "Nextclade_pango",
                    "title": "Nextclade Pango Lineage",
                    "type": "categorical"
                },
                {
                    "key": "S1_mutations",
                    "title": "S1 Mutations",
                    "type": "continuous"
                },
                {
                    "key": "immune_escape",
                    "title": "Immune Escape vs BA.2",
                    "type": "continuous"
                },
                {
                    "key": "ace2_binding",
                    "title": "ACE2 binding vs BA.2",
                    "type": "continuous"
                },
                {
                    "key": "mlr_lineage_fitness",
                    "title": "MLR lineage fitness",
                    "type": "continuous"
                },
                {
                    "key": "region",
                    "title": "Region",
                    "type": "categorical"
                },
                {
                    "key": "country",
                    "title": "Country",
                    "type": "categorical"
                },
                {
                    "key": "division",
                    "title": "Admin Division",
                    "type": "categorical"
                },
                location_coloring,
                {
                    "key": "host",
                    "title": "Host",
                    "type": "categorical"
                },
                {
                    "key": "author",
                    "title": "Authors",
                    "type": "categorical"
                },
                {
                    "key": "originating_lab",
                    "title": "Originating Lab",
                    "type": "categorical"
                },
                {
                    "key": "submitting_lab",
                    "title": "Submitting Lab",
                    "type": "categorical"
                },
                {
                    "key": "recency",
                    "title": "Submission Date",
                    "type": "categorical"
                },
                {
                    "key": "epiweek",
                    "title": "Epiweek (CDC)",
                    "type": "categorical"
                },
                gisaid_epi_isl_coloring,
                {
                    "key": "genbank_accession",
                    "type": "categorical"
                }
            ],
            "geo_resolutions": [
                "region",
                "country",
                "division"
            ],
            "display_defaults": {
                "color_by": "clade_membership",
                "distance_measure": "num_date",
                "geo_resolution": default_geo_resolution,
                "map_triplicate": default_map_triplicate,
                "branch_label": "clade",
                "transmission_lines": False
            },
            "filters": [
                "clade_membership",
                "emerging_lineage",
                "pango_lineage",
                "Nextclade_pango",
                "region",
                "country",
                "division",
                location_filter,
                "host",
                "author",
                originating_lab_filter,
                submitting_lab_filter,
                "recency",
            ],
            "panels": [
                "tree",
                "map",
                "entropy",
                "frequencies"
            ]
        }

        ## Prune out None values
        data['colorings'] = [c for c in data['colorings'] if c!=None]
        data['filters'] = [f for f in data['filters'] if f!=None]

        with open(output[0], 'w') as fh:
            json.dump(data, fh, indent=2)

#
# Rules for custom auspice exports for the Nextstrain team.
#

rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.finalize.output.auspice_json,
        tip_frequencies_json = rules.finalize.output.tip_frequency_json,
        root_sequence_json = rules.finalize.output.root_sequence_json
    output:
        dated_auspice_json = "auspice/{prefix}_{build_name}_{date}.json",
        dated_tip_frequencies_json = "auspice/{prefix}_{build_name}_{date}_tip-frequencies.json",
        dated_root_sequence_json = "auspice/{prefix}_{build_name}_{date}_root-sequence.json"
    benchmark:
        "benchmarks/dated_json_{prefix}_{build_name}_{date}.txt"
    wildcard_constraints:
        # Allow build names to contain alphanumeric characters, underscores, and
        # hyphens but not special strings used for Nextstrain builds. Include
        # the user-defined prefix as a constraint, so Snakemake does not parse
        # parts of the actual build names as part of the prefix.
        prefix = re.escape(config["auspice_json_prefix"]),
        build_name = r'(?:[-a-zA-Z0-9_](?!tip-frequencies|root-sequence|\d{4}-\d{2}-\d{2}))+',
        date = r"\d{4}-\d{2}-\d{2}"
    conda: config["conda_environment"]
    shell:
        r"""
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        cp {input.root_sequence_json} {output.dated_root_sequence_json}
        """

#
# Deployment and error handlers, including Slack messaging integrations.
#

from os import environ

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


# We limit the number of concurrent `nextstrain deploy` jobs by using a custom
# resource, concurrent_deploys, for the rule deploy_single (below).  This is
# because each `nextstrain deploy` to production currently performs a
# CloudFront wildcard invalidation, and these invalidations have a concurrency
# limit of 15 per CloudFront distribution.¹
#
# Set a default for concurrent_deploys so workflow invocations don't have to
# specify it with --resources or profile config.  This is particularly
# important because using --resources (or the equivalent profile config)
# conflicts with `nextstrain build`'s automatic passing thru of allocated
# memory to Snakemake.²
#
# The default is a very conservative 2 because `nextstrain deploy` does not
# always wait for invalidation completion before exiting, and thus
# invalidations may still be running even though the Snakemake job that started
# it is not.  We also typically have multiple workflow runs occurring
# simultaneously (e.g. GISAID and Open), and the limit of 15 is effectively
# shared between them.
#   -trs, 1 Feb 2023
#
# ¹ <https://docs.aws.amazon.com/AmazonCloudFront/latest/DeveloperGuide/cloudfront-limits.html#limits-invalidations>
# ² <https://github.com/nextstrain/cli/blob/017c5380/nextstrain/cli/command/build.py#L149-L166>
workflow.global_resources.setdefault("concurrent_deploys", 2)

rule deploy_single:
    input:
        # Replication of all_regions but with build_name as wildcard
        auspice_json = expand("auspice/{prefix}_{{build_name}}.json", prefix=config["auspice_json_prefix"]),
        tip_frequencies_json = expand("auspice/{prefix}_{{build_name}}_tip-frequencies.json", prefix=config["auspice_json_prefix"]),
        root_sequence_json = expand("auspice/{prefix}_{{build_name}}_root-sequence.json", prefix=config["auspice_json_prefix"]),
        dated_auspice_json = expand("auspice/{prefix}_{{build_name}}_{date}.json", prefix=config["auspice_json_prefix"], date=config.get("build_date", get_todays_date())),
        dated_tip_frequencies_json = expand("auspice/{prefix}_{{build_name}}_{date}_tip-frequencies.json", prefix=config["auspice_json_prefix"], date=config.get("build_date", get_todays_date())),
        dated_root_sequence_json = expand("auspice/{prefix}_{{build_name}}_{date}_root-sequence.json", prefix=config["auspice_json_prefix"], date=config.get("build_date", get_todays_date()))
    output:
        temp(touch("auspice/{build_name}.deployed"))
    params:
        deploy_url = config["deploy_url"],
        prefix = config["auspice_json_prefix"]
    resources:
        concurrent_deploys = 1
    benchmark:
        "benchmarks/deploy_single_{build_name}.txt"
    run:
        shell("nextstrain deploy {params.deploy_url:q} {input:q}")
        if params.deploy_url == "s3://nextstrain-data":
            slack_message = f"Deployed {wildcards.build_name} to https://nextstrain.org/{(params.prefix + '_' + wildcards.build_name).replace('_', '/')}"
        else:
            slack_message = f"Deployed {wildcards.build_name} to {params.deploy_url}"
        send_slack_message(slack_message)

rule deploy:
    input:
        expand("auspice/{build_name}.deployed", build_name=BUILD_NAMES),
    params:
        slack_message = f"Deployed to {config['deploy_url']} {deploy_origin}",
    benchmark:
        "benchmarks/deploy.txt"
    run:
        send_slack_message(params.slack_message)

rule upload:
    message: "Uploading intermediate files for specified origins to {params.s3_bucket}"
    input:
        unpack(_get_upload_inputs)
    output:
        touch("results/upload.done")
    params:
        s3_bucket = config["S3_DST_BUCKET"],
    log:
        "logs/upload.txt"
    benchmark:
        "benchmarks/upload.txt"
    run:
        message = "The following files have been updated (unless they were identical):"
        for remote, local in input.items():
            shell("./scripts/upload-to-s3 {local:q} s3://{params.s3_bucket:q}/{remote:q} | tee -a {log:q}")
            message += f"\n\ts3://{params.s3_bucket}/{remote}"
        send_slack_message(message)

storage = PersistentDict("slack")

def send_slack_message(message, broadcast=False):
    """
    Sends slack messages notifying us of the pipeline's progress. Messages will be
    threaded, with the first message being the parent message. Important messages can be
    broadcast: they will be part of the thread but also sent to the channel.
    """
    ## Note: this cannot currently send files, but this would be easy to add.
    ## Slack docs: https://api.slack.com/methods/files.upload

    if not config.get("slack_token", None) or not config.get("slack_channel", None):
        print("Cannot send slack message as the config does not define a channel and/or token.")
        return

    headers = {
        'Content-type': 'application/json',
        'authorization': f"Bearer {config['slack_token']}",
        'Accept': 'application/vnd.github.v3+json'
    }
    data = {
        "channel": config["slack_channel"],
        "text": message,
    }

    # if slack_thread_ts has been stored, then there is a parent message, so we thread this messaege
    try:
        data["thread_ts"]=str(storage.fetch("slack_thread_ts"))
        if broadcast:
            data["reply_broadcast"]=True
    except NoSuchEntryError:
        pass

    response = requests.post("https://slack.com/api/chat.postMessage", headers=headers, data=json.dumps(data))
    response.raise_for_status()
    storage.store_if_not_present("slack_thread_ts", response.json()["ts"])

# onstart handler will be executed before the workflow starts.
onstart:
    message = [
        "Pipeline starting which will run the phylogenetics and upload build assets",
        f"Deployed {deploy_origin}"
    ]
    if environ.get("AWS_BATCH_JOB_ID"):
        message.append(f" (<https://console.aws.amazon.com/batch/v2/home?region=us-east-1#jobs/detail/{environ['AWS_BATCH_JOB_ID']}|AWS console link>)")
    message.append("Further results will appear in this 🧵.")
    send_slack_message(". ".join(message))

# onsuccess handler is executed if the workflow finished without error.
onsuccess:
    message = "✅ This pipeline has successfully finished 🎉"
    send_slack_message(message)
    storage.clear() # clear any persistent storage

# onerror handler is executed if the workflow finished with an error.
onerror:
    message = "❌ This pipeline has FAILED 😞. Please see linked thread for more information."
    send_slack_message(message, broadcast=True)
    storage.clear() # clear any persistent storage
