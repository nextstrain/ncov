"""
calculate priorties from index and proximities
"""
import argparse
from random import shuffle
from collections import defaultdict
import numpy as np
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequence-index", type=str, required=True, help="sequence index file")
    parser.add_argument("--proximities", type = str, required=True, help="tsv file with proximities")
    parser.add_argument("--Nweight", type = float, default=0.003, required=False, help="parameterizes de-prioritization of incomplete sequences")
    parser.add_argument("--crowding-penalty", type = float, default=0.01, required=False, help="parameterizes how priorities decrease when there is many very similar sequences")
    parser.add_argument("--output", type=str, required=True, help="tsv file with the priorities")
    args = parser.parse_args()

    proximities = pd.read_csv(args.proximities, sep='\t', index_col=0)
    index = pd.read_csv(args.sequence_index, sep='\t', index_col=0)
    combined = pd.concat([proximities, index], axis=1)

    closest_matches = combined.groupby('closest strain')
    candidates = {}
    for focal_seq, seqs in closest_matches.groups.items():
        tmp = combined.loc[seqs, ["distance", "N"]]
        # penalize larger distances and more undetermined sites. 1/args.Nweight are 'as bad' as one extra mutation
        tmp["priority"] = -tmp.distance - tmp.N*args.Nweight
        name_prior = [(name, d.priority) for name, d in tmp.iterrows()]
        shuffle(name_prior)
        candidates[focal_seq] = sorted(name_prior, key=lambda x:x[1])

    # export priorities
    crowding = args.crowding_penalty
    with open(args.output, 'w') as fh:
        # loop over lists of sequences that are closest to particular focal sequences
        for cs in candidates.values():
            # these sets have been shuffled -- reduce priorities in this shuffled random order
            for i, (name, pr) in enumerate(cs):
                fh.write(f"{name}\t{pr-i*crowding:1.2f}\n")
