from augur.io import read_metadata
from augur.utils import write_json
import requests
import json
import pandas as pd
import argparse
import math

# Set up argument parser
parser = argparse.ArgumentParser(description="Process metadata and growth advantage data.")
parser.add_argument("--metadata", required=True, help="Path to the metadata file (TSV or compressed .tsv.xz format).")
parser.add_argument("--metadata-id-columns", default=["strain", "name", "Virus name"], nargs="+", help="List of columns to use as identifiers in the metadata file.")
parser.add_argument("--metadata-clade-attribute", default="Nextclade_pango", help="Matched attribute to MLR variants.")
parser.add_argument("--mlr-url", default="https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/global/mlr/latest_results.json", help="URL to fetch the forecasts JSON data.")
parser.add_argument("--output-node-data", required=True, help="Path to save the output JSON node data.")

args = parser.parse_args()

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
        print(f"Error fetching the JSON file: {e}")
        return None

try:
    # Fetch the growth advantages
    growth_advantages = fetch_growth_advantages(args.mlr_url)

    # Load the local metadata file
    metadata_file = args.metadata
    metadata = read_metadata(
        metadata_file,
        id_columns=args.metadata_id_columns
    )

    # Match Nextclade_pango entries to the growth advantage
    if growth_advantages:
        metadata[args.metadata_clade_attribute] = metadata[args.metadata_clade_attribute].map(growth_advantages)
    else:
        metadata[args.metadata_clade_attribute] = math.nan

    # Output rows with matched data
    print(metadata.head())  # Display the first few rows as an example

    # Create a node data object with growth advantages
    node_data = {}
    for index, record in metadata.iterrows():
        node_data[index] = {
            "mlr_lineage_fitness": record[args.metadata_clade_attribute]
        }

    # Save node data
    write_json({"nodes": node_data}, args.output_node_data)

except FileNotFoundError as e:
    print(f"Error reading metadata file: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
