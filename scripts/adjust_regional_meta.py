"""
Very simple script just to change 'division' & 'location' to 'country' for regions outwith focal region
"""

import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--region", type=str, required=True, help="region for which to 'mask' at division level & below")
    parser.add_argument("--output", type=str, required=True, help="new metadata")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, delimiter='\t')
    metadata.insert(11, 'focal_country', metadata['country'])
    metadata.loc[metadata.region != args.region, 'division'] = ""
    metadata.loc[metadata.region != args.region, 'location'] = ""
    metadata.loc[metadata.region != args.region, 'focal_country'] = ""

    metadata.to_csv(args.output, index=False, sep="\t")
