#!/usr/bin/env python3
import argparse
from augur.utils import write_json
import epiweeks
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Calculate epiweeks for dates in the given metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="metadata with a 'date' column")
    parser.add_argument("--attribute-name", default="epiweek", help="name to store annotations of epiweeks in JSON output")
    parser.add_argument("--output-node-data", required=True, help="node data JSON with epiweek annotations")

    args = parser.parse_args()

    # Read metadata with pandas because Augur's read_metadata utility does not
    # support metadata without a "strain" or "name" field.
    metadata = pd.read_csv(
        args.metadata,
        sep=None,
        engine="python",
        skipinitialspace=True,
        dtype={
            "strain": "string",
            "name": "string",
        }
    ).fillna("")

    # Find records with unambiguous dates.
    metadata_with_dates = metadata.loc[~metadata["date"].str.contains("X"), ["strain", "date"]].copy()

    # Convert date strings to timestamps.
    metadata_with_dates["date"] = pd.to_datetime(metadata_with_dates["date"])

    # Calculate epiweeks from date objects as a new annotation.
    metadata_with_dates["epiweek"] = metadata_with_dates["date"].apply(lambda date: epiweeks.Week.fromdate(date).cdcformat())

    # Create a node data object with epiweeks.
    node_data = {}
    for record in metadata_with_dates.to_dict(orient="records"):
        node_data[record["strain"]] = {
            args.attribute_name: record["epiweek"],
        }

    # Save node data.
    write_json({"nodes": node_data}, args.output_node_data)
