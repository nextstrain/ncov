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

import requests
import json
from workflow.lib.persistent_dict import PersistentDict, NoSuchEntryError

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date

rule all_regions:
    input:
        auspice_json = expand("auspice/{prefix}_{build_name}.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES),
        tip_frequencies_json = expand("auspice/{prefix}_{build_name}_tip-frequencies.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES),
        dated_auspice_json = expand("auspice/{prefix}_{build_name}_{date}.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES, date=config.get("build_date", get_todays_date())),
        dated_tip_frequencies_json = expand("auspice/{prefix}_{build_name}_{date}_tip-frequencies.json", prefix=config["auspice_json_prefix"], build_name=BUILD_NAMES, date=config.get("build_date", get_todays_date()))

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
        mem_mb=lambda wildcards, input: 5 * int(input.metadata.size / 1024 / 1024)
    conda: config["conda_environment"]
    shell:
        """
        python3 ./scripts/check_missing_locations.py \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --latlong {input.lat_longs}
        """


rule mutation_summary:
    message: "Summarizing {input.alignment}"
    input:
        alignment = rules.align.output.alignment,
        insertions = rules.align.output.insertions,
        translations = rules.align.output.translations,
        reference = config["files"]["alignment_reference"],
        genemap = config["files"]["annotation"]
    output:
        mutation_summary = "results/mutation_summary_{origin}.tsv.xz"
    log:
        "logs/mutation_summary_{origin}.txt"
    benchmark:
        "benchmarks/mutation_summary_{origin}.txt"
    params:
        outdir = "results/translations",
        basename = "seqs_{origin}",
        genes=config["genes"],
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mutation_summary.py \
            --alignment {input.alignment} \
            --insertions {input.insertions} \
            --directory {params.outdir} \
            --basename {params.basename} \
            --reference {input.reference} \
            --genes {params.genes:q} \
            --genemap {input.genemap} \
            --output {output.mutation_summary} 2>&1 | tee {log}
        """


#
# Rules for custom auspice exports for the Nextstrain team.
#

rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.finalize.output.auspice_json,
        tip_frequencies_json = rules.include_hcov19_prefix.output.tip_frequencies
    output:
        dated_auspice_json = "auspice/{prefix}_{build_name}_{date}.json",
        dated_tip_frequencies_json = "auspice/{prefix}_{build_name}_{date}_tip-frequencies.json"
    benchmark:
        "benchmarks/dated_json_{prefix}_{build_name}_{date}.txt"
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

rule deploy:
    input:
        *rules.all_regions.input
    params:
        slack_message = f"Deployed to {config['deploy_url']} {deploy_origin}",
        deploy_url = config["deploy_url"]
    benchmark:
        "benchmarks/deploy.txt"
    run:
        shell("nextstrain deploy {params.deploy_url:q} {input:q}")
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
