import argparse, json
from Bio import Phylo, SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
from augur.distance import read_distance_map

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--map', nargs='+')
    parser.add_argument('--attribute-name', nargs='+')
    parser.add_argument('--simple-map', action='store_true')
    args = parser.parse_args()

    if type(args.map)==str:
        args.map = [args.map]

    if type(args.attribute_name)==str:
        args.attribute_name = [args.attribute_name]

    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0)['S'].fillna('')
    positions = {}
    positions_with_muts = {}
    for s,v in mutations.items():
        positions[s] = [int(x[1:-1]) for x in v.split(',')] if v else []
        positions_with_muts[s] = [(int(x[1:-1]), x[0], x[-1]) for x in v.split(',')] if v else []
        if len(positions[s])>30:
            print('excluding',s)
            positions[s]=[]

    L=1400
    scores = defaultdict(dict)

    if args.simple_map:
        for map_file, attribute in zip(args.map, args.attribute_name):
            with open(map_file) as fh:
                map = json.load(fh)
            weights = map['default']*np.ones(L)
            for k,v in map["map"]['S'].items():
                weights[int(k)-1] = v

            for s,pos in positions.items():
                scores[s][attribute] = np.sum(weights[pos])
    else:
        for map_file, attribute in zip(args.map, args.attribute_name):
            dmap = read_distance_map(map_file)
            for s, muts in positions_with_muts.items():
                val = 0
                for pos, wt, m in muts:
                    if pos in dmap["map"]:
                        val += dmap["map"][pos].get((wt,m),0)
                scores[s][attribute] = val

    d = pd.DataFrame(scores).T

    all_scores = pd.concat([mutations, d], axis=1)
    all_scores["sum"] = d.sum(axis=1)
    all_scores["number"] = (d>0).sum(axis=1)


    thres = 20
    for mut_count in [2,3]:
        mult_muts = all_scores.number==mut_count
        mult_hits = all_scores.loc[mult_muts,:]
        print(f"\n# Substantial escape in {mut_count} class vaccines")
        for r, row in mult_hits.iterrows():
            if row['sum']>thres:
                print(f"{r}\t{row.S}\t{row.c1}\t{row.c2}\t{row.c3}\t{row['sum']}")

