import argparse, json
from Bio import Phylo, SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
from augur.distance import read_distance_map
from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--min-date', type=str)
    parser.add_argument('--min-sub-date', type=str)
    parser.add_argument('--country')
    parser.add_argument('--region')
    args = parser.parse_args()

    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0)['S'].fillna('')
    meta = pd.read_csv(args.metadata, sep='\t', index_col=0).fillna('')
    data = pd.concat([mutations, meta.loc[:, ["date", "date_submitted", "country", "region"]]], axis=1).loc[mutations.index]

    if args.min_date:
        data = data.loc[data.date>args.min_date]

    if args.min_sub_date:
        data = data.loc[data.date_submitted>args.min_sub_date]

    if args.country:
        data = data.loc[data.country==args.country]

    if args.region:
        data = data.loc[data.region==args.region]
        

    positions = {}
    positions_with_muts = {}
    for s,v in data.iterrows():
        positions[s] = [int(x[1:-1]) for x in v['S'].split(',')] if v['S'] else []
        positions_with_muts[s] = [(int(x[1:-1]), x[0], x[-1]) for x in v['S'].split(',')] if v['S'] else []
        if len(positions[s])>30:
            print('excluding',s)
            positions[s]=[]
            positions_with_muts[s]=[]

    combos = {}
    for combi in powerset([18,417,484,501,681]):
        if len(combi)==0:
            continue
        tmp = []
        for k,v in positions:
            if all([x in v for x in combi]):
                tmp.append(k)
        combos[combi] = tmp


