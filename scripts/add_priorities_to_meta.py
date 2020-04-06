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
        description="Add columns for priorities of sequences relative to diff focal regions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--priorities", type=str, nargs="+", required=True, help="priorities files")
    parser.add_argument("--config", type=str, help="config file to modify")
    parser.add_argument("--output-meta", type=str, required=True, help="adjusted metadata")
    parser.add_argument("--output-config", type=str, help="modified config")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')
    with open(args.config) as fh:
        input_json = json.load(fh)

    for priority_file in args.priorities:
        p_f = priority_file.replace(".tsv", "")
        region = p_f.split("_")[2]
        column_name = "".join(["priorities_",region])
        
        with open(priority_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            priors = {r[0]: r[1] for r in reader if len(r)>1}

        assign_priors =  [priors[st] if st in priors else "" for st in metadata.strain]

        metadata.insert(11, column_name, assign_priors)
        input_json['colorings'].append({'key': column_name, 'type': 'continuous'})

    metadata.to_csv(args.output_meta, index=False, sep="\t")

    with open(args.output_config, 'w') as fh:
        json.dump(input_json, fh, indent=2)
    