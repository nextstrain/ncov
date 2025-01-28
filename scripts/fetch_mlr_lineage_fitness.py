import argparse
import json
import pandas as pd
import requests
import sys
from augur.io import read_metadata
from augur.utils import write_json

# This script is currently assuming a match on lineage fitness this uses
# https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/global/mlr/latest_results.json
# that backs the live estimates on https://nextstrain.org/sars-cov-2/forecasts
# This uses "Nextclade_pango" metadata label to derive sequence counts from
# GISAID data and estimate relative growth advantages across collapsed Pango
# lineages. It will be most relevant for 1m, 2m and 6m builds, but is not at all
# broken for the all-time builds. It would be possible to swap this to key on
# clade instead, but I think the greater detail of lineages is better in this case

def fetch_growth_advantages(mlr_url):
    try:
        response = requests.get(mlr_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        json_data = response.json()  # Parse the JSON content
        data = json_data["data"]

        growth_advantages = {}
        for entry in data:
            if all(key in entry for key in ["location", "site", "variant", "value", "ps"]):
                if entry["location"] == "hierarchical" and entry["site"] == "ga" and entry["ps"] == "median":
                    growth_advantages[entry["variant"]] = entry["value"]
        return growth_advantages
    except Exception as e:
        print(f"Error fetching the JSON file: {e}", file=sys.stderr)
        return None

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Fetch MLR lineage fitness and match to strain-level metadata")
    parser.add_argument("--metadata", required=True, help="Path to the metadata TSV")
    parser.add_argument("--metadata-id-columns", default=["strain", "name", "Virus name"], nargs="+", help="List of columns to use as identifiers in the metadata file")
    parser.add_argument("--metadata-clade-attribute", default="Nextclade_pango", help="Matched attribute to MLR variants")
    parser.add_argument("--mlr-url", default="https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/global/mlr/latest_results.json", help="URL to fetch the forecasts JSON data.")
    parser.add_argument("--output-node-data", required=True, help="Path to save the output JSON node data.")

    args = parser.parse_args()

    # Fetch the remote growth advantages
    growth_advantages = fetch_growth_advantages(args.mlr_url)

    # Load the local metadata
    try:
        metadata = read_metadata(
            args.metadata,
            id_columns=args.metadata_id_columns
        )
    except FileNotFoundError as e:
        print(f"Error reading metadata file: {e}", file=sys.stderr)

    # Match Nextclade_pango entries in metadata to the fetched growth advantages
    if growth_advantages:
        metadata[args.metadata_clade_attribute] = metadata[args.metadata_clade_attribute].map(growth_advantages)
    else:
        metadata[args.metadata_clade_attribute] = None

    # Create a node data object with growth advantages
    node_data = {}
    for index, record in metadata.iterrows():
        node_data[index] = {
            "mlr_lineage_fitness": record[args.metadata_clade_attribute] if pd.notna(record[args.metadata_clade_attribute]) else None
        }

    # Save node data
    write_json({"nodes": node_data}, args.output_node_data)

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
