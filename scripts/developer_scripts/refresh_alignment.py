import argparse
import boto3
import datetime
import os
import subprocess
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Trigger AWS Batch job to refresh alignment if stale",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--source', type=str, default="gisaid", help="source to compare, gisaid vs open")
    args = parser.parse_args()

    # confirm that script is run from correct directory
    current_directory = str(os.getcwd()).split("/")[-1]
    if (current_directory != "ncov"):
        print("Script must be run from the parent ncov/ directory")
        sys.exit()

    # confirm that Nextstrain CLI is available
    try:
        subprocess.call(["nextstrain"], stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        print("Nextstrain CLI must be available")
        sys.exit()

    s3 = boto3.resource('s3')
    bucket = 'nextstrain-ncov-private'

    sequences_date = s3.ObjectSummary(bucket, 'sequences.fasta.xz').last_modified
    print("sequences_date", sequences_date)

    filtered_date = s3.ObjectSummary(bucket, 'filtered.fasta.xz').last_modified
    print("filtered_date", filtered_date)

    difference = sequences_date - filtered_date
    hours_stale = divmod(difference.total_seconds(), 3600)[0]
    print("hours_stale", hours_stale)

    # trigger AWS Batch job if more stale than 6 hours
    # requesting results/filtered_gisaid.fasta.xz will also generate:
    # - results/aligned_gisaid.fasta.xz
    # - results/masked_gisaid.fasta.xz
    if (hours_stale > 1):
        call = ['nextstrain', 'build', '--aws-batch', '--cpus', '8', '--memory', '14GiB', '--detach', \
            '.', '--jobs', '8', '--profile', 'nextstrain_profiles/nextstrain-gisaid', 'results/filtered_gisaid.fasta.xz']
        print(' '.join(call))
