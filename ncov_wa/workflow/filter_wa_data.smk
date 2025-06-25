rule pull_full_data:
    output:
        full_metadata="data/full_metadata.tsv.gz",
        full_sequences="data/full_sequences.fasta.xz"
    params:
        full_metadata_url="https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz",
        full_sequences_url="https://data.nextstrain.org/files/ncov/open/sequences.fasta.xz"
    shell:
        """
        bash ncov_wa/scripts/pull_full_data.sh {params.full_metadata_url} {params.full_sequences_url} {output.full_metadata} {output.full_sequences}
        """


rule filter_wa_sequences_step:
    input:
        fasta="data/full_sequences.fasta.xz"
    output:
        filtered_wa_sequences="data/filtered_wa_sequences.fasta"
    shell:
        """
        bash ncov_wa/scripts/filter_wa_sequences.sh {input.fasta} {output.filtered_wa_sequences}
        """


rule filter_wa_metadata_step:
    input:
        filtered_wa_sequences="data/filtered_wa_sequences.fasta",
        metadata="data/full_metadata.tsv.gz",
        headers="ncov_wa/data/headers.tsv"
    output:
        temp="ncov_wa/data/tmp.tsv",
        filtered_wa_metadata="data/filtered_wa_metadata.tsv"
    shell:
        """
        bash ncov_wa/scripts/filter_wa_metadata.sh {input.filtered_wa_sequences} {input.metadata} {input.headers} {output.temp} {output.filtered_wa_metadata}
        """

rule add_county_metadata:
    input:
        filtered_wa_metadata="data/filtered_wa_metadata.tsv",
        county_metadata="ncov_wa/data/county_metadata.csv"
    output:
        wa_metadata="data/wa-metadata.tsv"
    shell:
        """
        python3 ncov_wa/scripts/wa-nextstrain-update-location-genbank.py {input.filtered_wa_metadata} {input.county_metadata} {output.wa_metadata}
        """


rule compress_files:
    input:
        filtered_wa_sequences="data/filtered_wa_sequences.fasta",
        wa_metadata="data/wa-metadata.tsv"
    output:
        compressed_output="data/wa-sequences.tar.xz"
    shell:
        """
        tar -cJf {output.compressed_output} {input.filtered_wa_sequences} {input.wa_metadata}
        """

# to set build to use WA data from the past six months
from dateutil import relativedelta

# Calculate dates
d = date.today()
six_m = d - relativedelta.relativedelta(months=6)

# Set earliest_date & latest_date in builds
if "ncov_wa_six_mon" in config["builds"]:
    config["builds"]["ncov_wa_six_mon"]["earliest_date"]= six_m.strftime('%Y-%m-%d')


ruleorder:  compress_files > filter_wa_metadata_step > filter_wa_sequences_step > pull_full_data
