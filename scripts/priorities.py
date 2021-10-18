"""
calculate priorties from index and proximities
"""
import argparse
from random import shuffle
from collections import defaultdict
import numpy as np
import pandas as pd

def create_priorities(sequence_index_path, proximities_path, output_path, Nweight=0.003, crowding_penalty=0.01):
    """
    calculate priorties from index and proximities
    Parameters
    ----------
    sequence_index_path : string
        Path to sequence index file
    proximities_path : string
        path to tsv file with proximities
    output_path : string
        path to TSV file with the priorities
    Nweight : float, default: 0.003
        parameterizes de-prioritization of incomplete sequences
    crowding_penalty : float, default: 0.01
        parameterizes how priorities decrease when there is many very similar sequences
    Returns
    -------
    None
    """

    proximities = pd.read_csv(proximities_path, sep='\t', index_col=0)
    index = pd.read_csv(sequence_index_path, sep='\t', index_col=0)
    combined = pd.concat([proximities, index], axis=1)

    closest_matches = combined.groupby('closest strain')
    candidates = {}
    for focal_seq, seqs in closest_matches.groups.items():
        tmp = combined.loc[seqs, ["distance", "N"]]
        # penalize larger distances and more undetermined sites. 1/args.Nweight are 'as bad' as one extra mutation
        tmp["priority"] = -tmp.distance - tmp.N*Nweight
        name_prior = [(name, d.priority) for name, d in tmp.iterrows()]
        shuffle(name_prior)
        candidates[focal_seq] = sorted(name_prior, key=lambda x:x[1])

    # export priorities
    crowding = crowding_penalty
    with open(output_path, 'w') as fh:
        # loop over lists of sequences that are closest to particular focal sequences
        for cs in candidates.values():
            # these sets have been shuffled -- reduce priorities in this shuffled random order
            for i, (name, pr) in enumerate(cs):
                fh.write(f"{name}\t{pr-i*crowding:1.2f}\n")


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
    create_priorities(
        args.sequence_index,
        args.proximities,
        args.output,
        args.Nweight,
        args.crowding_penalty
    )
