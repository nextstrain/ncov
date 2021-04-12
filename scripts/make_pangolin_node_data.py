"""
Translate pangolineages from CSV -> JSON for node_data
Note: this should arguably live instead as part of `combine_metadata`,
but this gets particularly complex given the new multiple-inputs logic.
So, for now, following the initial suggestion in the issue.
"""

import argparse
import pandas as pd
import csv
import json
from augur.utils import write_json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create node data for assigned pangolin lineages",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pangolineages", type=str, required=True, help="pangolineages.csv")
    parser.add_argument("--node_data_outfile", type=str, help="pangolineages.json")
    parser.add_argument("--attribute_name", default="pango_lineage_local", help="attribute name for pangolin lineage annotations in the output JSON")
    args = parser.parse_args()

    pangolineages = pd.read_csv(args.pangolineages)

    node_data = {
    "nodes": {
    row['taxon']: {args.attribute_name: row['lineage']} for idx, row in pangolineages.iterrows()
        }
    }

    write_json(node_data, args.node_data_outfile)
