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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create node data for assigned pangolin lineages",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pangolineages", type=str, required=True, help="pangolineages.csv")
    parser.add_argument("--node_data_outfile", type=str, help="pangolineages.json")
    args = parser.parse_args()
    print('INPUT TO MAKE NODE DATA', '\n\n', args.pangolineages, '\n\n', args.node_data_outfile)

    pangolineages = pd.read_csv(args.pangolineages)
    print(pangolineages.head())
    node_data = {
    "nodes": {
    row['taxon']: row['lineage'] for idx, row in pangolineages.iterrows()
        }
    }

    # input_json['colorings'].append({'key': 'pangolin lineage', 'type': 'categorical'})

    with open(args.node_data_outfile, 'w') as fh:
        json.dump(node_data, fh, indent=2)
