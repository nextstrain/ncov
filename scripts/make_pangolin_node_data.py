"""
Add column to metadata with the priorities of 'context' sequences
relative to the 'focal' samples
"""

import argparse
import pandas as pd
import csv
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create node data for assigned pangolin lineages",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pangolineages", type=str, nargs="1", required=True, help="pangolineages.csv")
    parser.add_argument("--node_data_outfile", type=str, help="pangolineages.json")
    args = parser.parse_args()

    pangolineages = pd.read_csv(args.pangolineages, sep='\t')
    node_data = {
    "nodes": {
    row['Taxon']: row['Lineage'] for idx, row in pangolineages.iterrows()
        }
    }

    # input_json['colorings'].append({'key': 'pangolin lineage', 'type': 'categorical'})

    with open(args.node_data_outfile, 'w') as fh:
        json.dump(node_data, fh, indent=2)
