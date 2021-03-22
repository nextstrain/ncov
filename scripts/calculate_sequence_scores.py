import argparse, json
from Bio import Phylo, SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--map', nargs='+')
    parser.add_argument('--attribute-name', nargs='+')
    args = parser.parse_args()

    if type(args.map)==str:
        args.map = [args.map]

    if type(args.attribute_name)==str:
        args.attribute_name = [args.attribute_name]

    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0)['S'].fillna('')
    positions = {}
    for s,v in mutations.items():
        positions[s] = [int(x[1:-1])-1 for x in v.split(',')] if v else []
        if len(positions[s])>30:
            print('excluding',s)
            positions[s]=[]

    L=1400
    scores = defaultdict(dict)

    for map_file, attribute in zip(args.map, args.attribute_name):
        with open(map_file) as fh:
            map = json.load(fh)
        weights = map['default']*np.ones(L)
        for k,v in map["map"]['S'].items():
            weights[int(k)-1] = v

        for s,pos in positions.items():
            scores[s][attribute] = np.sum(weights[pos])

    d = pd.DataFrame(scores).T

    all_scores = pd.concat([mutations, d], axis=1)
    all_scores["sum"] = d.sum(axis=1)
    all_scores["number"] = (d>0).sum(axis=1)

    triple_muts = all_scores.number==3
    all_scores.loc[triple_muts,:]

    # with open(args.output, 'w') as fh:
    #     json.dump({"nodes":node_data}, fh)
