#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
from augur.utils import write_json
import epiweeks
import pandas as pd
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Calculate epiweeks for dates in the given metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="metadata with a 'date' column")
    parser.add_argument("--metadata-id-columns", default=["strain", "name", "Virus name"], nargs="+", help="names of valid metadata columns containing identifier information like 'strain' or 'name'")
    parser.add_argument("--attribute-name", default="epiweek", help="name to store annotations of epiweeks in JSON output")
    parser.add_argument("--output-node-data", required=True, help="node data JSON with epiweek annotations")

    args = parser.parse_args()

    metadata = read_metadata(
        args.metadata,
        id_columns=args.metadata_id_columns,
    )

    # Find records with unambiguous dates. These must be complete date-like
    # records in YYYY-MM-DD format.
    date_pattern = re.compile(r"^\d{4}-\d{2}-\d{2}$")
    has_complete_date = metadata["date"].astype(str).apply(lambda date: date_pattern.match(date) is not None)
    metadata_with_dates = metadata.loc[has_complete_date, ["date"]].copy()

    # Convert date strings to timestamps.
    metadata_with_dates["date"] = pd.to_datetime(metadata_with_dates["date"])

    # Calculate epiweeks from date objects as a new annotation.
    metadata_with_dates["epiweek"] = metadata_with_dates["date"].apply(lambda date: epiweeks.Week.fromdate(date).cdcformat())

    # Create a node data object with epiweeks.
    node_data = {}
    for index, record in metadata_with_dates.iterrows():
        node_data[index] = {
            args.attribute_name: record["epiweek"],
        }

    # Save node data.
    write_json({"nodes": node_data}, args.output_node_data)
