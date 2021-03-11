"""
Run a number QC checks on an alignment. these involve divergence from a reference, clusters of mutations, and completeness
"""
import argparse
from collections import defaultdict
import numpy as np
from augur.utils import read_metadata
from datetime import datetime

tmrca = datetime(2019, 12, 1).toordinal()

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--output-exclusion-list", type=str, required=True, help="Output to-be-reviewed addition to exclude.txt")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')

    clock_deviation = np.array([float(x) if isfloat(x) else np.nan for x in metadata.clock_deviation])

    to_exclude = np.abs(clock_deviation)>15 | (metadata.snp_clusters>1)

    # write out file with sequences flagged for exclusion sorted by date
    with open(args.output_exclusion_list, 'w') as excl:
        for s in metadata.loc[to_exclude,'strain']:
            excl.write(f'{s}\n')

