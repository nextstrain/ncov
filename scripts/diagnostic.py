"""
Run a number QC checks on an alignment. these involve divergence from a reference, clusters of mutations, and completeness
"""
import argparse
import numpy as np
import pandas as pd

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
    parser.add_argument("--clock-filter", type=float, default=15, help="maximal allowed deviation from the molecular clock")
    parser.add_argument("--snp-clusters", type=int, default=1, help="maximal allowed SNP clusters (as defined by nextclade)")
    parser.add_argument("--output-exclusion-list", type=str, required=True, help="Output to-be-reviewed addition to exclude.txt")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')

    if "clock_deviation" in metadata.columns:
        clock_deviation = np.array([float(x) if isfloat(x) else np.nan for x in metadata.clock_deviation])
    else:
        clock_deviation = np.zeros(len(metadata), dtype=bool)

    if "snp_clusters" in metadata.columns:
        snp_clusters = np.array([float(x) if isfloat(x) else np.nan for x in metadata.snp_clusters])
    else:
        snp_clusters = np.zeros(len(metadata), dtype=bool)
        
    to_exclude = (np.abs(clock_deviation)>args.clock_filter) | (snp_clusters>args.snp_clusters)

    # write out file with sequences flagged for exclusion sorted by date
    with open(args.output_exclusion_list, 'w') as excl:
        for s in metadata.loc[to_exclude,'strain']:
            excl.write(f'{s}\n')

