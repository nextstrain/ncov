import argparse, json
from Bio import Phylo, SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
from augur.distance import read_distance_map
from binding_calculator import BindingCalculator

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--mutation-summary', type=str)
    parser.add_argument('--metadata', type=str)
    parser.add_argument('--map', nargs='+')
    parser.add_argument('--attribute-name', nargs='+')
    parser.add_argument('--min-date', type=str)
    parser.add_argument('--min-sub-date', type=str)
    parser.add_argument('--country')
    parser.add_argument('--region')
    parser.add_argument('--simple-map', action='store_true')
    parser.add_argument('--threshold', type=float, default=20)
    args = parser.parse_args()

    if type(args.map)==str:
        args.map = [args.map]

    if type(args.attribute_name)==str:
        args.attribute_name = [args.attribute_name]

    bc = BindingCalculator()

    L=700
    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0)['S'].fillna('')
    mutations = mutations.apply(lambda x:','.join([y for y in x.split(',') if y and int(y[1:-1])<L]))
    meta = pd.read_csv(args.metadata, sep='\t', index_col=0).fillna('')
    all_scores = pd.concat([mutations, meta.loc[:, ["date", "date_submitted", "country", "region", "Nextstrain_clade", "pango_lineage"]]], axis=1).loc[mutations.index]

    if args.min_date:
        all_scores = all_scores.loc[all_scores.date>args.min_date]

    if args.min_sub_date:
        all_scores = all_scores.loc[all_scores.date_submitted>args.min_sub_date]

    if args.country:
        all_scores = all_scores.loc[all_scores.country==args.country]

    if args.region:
        all_scores = all_scores.loc[all_scores.region==args.region]
        
    positions = {}
    positions_with_muts = {}
    for s,v in all_scores.iterrows():
        positions[s] = [int(x[1:-1])-1 for x in v['S'].split(',')] if v['S'] else []
        positions_with_muts[s] = [(int(x[1:-1])-1, x[0], x[-1]) for x in v['S'].split(',')] if v['S'] else []
        if len(positions[s])>30:
            print('excluding',s)
            positions[s]=[]
            positions_with_muts[s]=[]

    scores = defaultdict(dict)

    if args.simple_map:
        for map_file, attribute in zip(args.map, args.attribute_name):
            with open(map_file) as fh:
                map = json.load(fh)
            weights = map['default']*np.ones(L)
            for k,v in map["map"]['S'].items():
                weights[int(k)-1] = v

            for s,pos in positions.items():
                scores[s][attribute] = np.sum(weights[[p for p in pos if p<L]])
    else:
        for map_file, attribute in zip(args.map, args.attribute_name):
            dmap = read_distance_map(map_file)
            for s, muts in positions_with_muts.items():
                val = 0
                for pos, wt, m in muts:
                    if pos in dmap["map"]['S']:
                        val += dmap["map"]['S'][pos].get((wt,m),0)
                scores[s][attribute] = val
    counter =  0
    for s,pos in positions.items():
        counter+=1
        if counter%1000==0:
            print(counter)
        scores[s]['binding_lost'] = 1-bc.binding_retained([p+1 for p in pos if p>330 and p<531])

    d = pd.DataFrame(scores).T
    all_scores = pd.concat([all_scores, d], axis=1)
    all_scores["sum"] = d.sum(axis=1)
    all_scores["number"] = (d>0).sum(axis=1)

    by_spike_groups = all_scores.groupby(by='S')
    by_spike_scores = by_spike_groups.mean()

    # this only works with the discrete maps, not the DMS ones
    # thres = args.threshold
    # for mut_count in [2,3,4]:
    #     mult_muts = all_scores.number==mut_count
    #     mult_hits = all_scores.loc[mult_muts,:]
    #     print(f"\n# Substantial escape in {mut_count} class epitopes")
    #     for r, row in mult_hits.iterrows():
    #         if row['sum']>thres:
    #             print(f"{r}\t{row.Nextstrain_clade}\t{row.pango_lineage}\t{row.S}\t{row.c1}\t{row.c2}\t{row.c3}\t{row['sum']}")


    nMax = 100
    # print the nMax sequences with the highest score for each map and the sum
    for a in args.attribute_name+['sum', 'binding_lost']:
        print("\n\n### ATTRIBUTE ",a)
        print(f"#{'strain':<40}\t{a}\t{'sum'}\tSpike mutations")
        for muts, row in by_spike_scores.sort_values(by=a)[-nMax:].iterrows():
            G = all_scores.loc[by_spike_groups.groups[muts]]
            if len(G)<5:
                continue
            print(f"{muts:<80}\t{G.date.min()}\t{G.date.max()}\t{len(G)}\t{G.Nextstrain_clade.mode()[0]}\t{G.pango_lineage.mode()[0]}\t{row[a]:1.2f}")



        # for r, row in all_scores.sort_values(by=a)[-nMax:].iterrows():
        #     print(f"{r:<40}\t{row.Nextstrain_clade}\t{row.pango_lineage}\t{row[a]:1.2f}\t{row['sum']:1.2f}\t{row.S}")


#major_mutations = mutations.apply(lambda x:','.join([y for y in x.split(',') if y and int(y[1:-1]) in [144,152,154,346,378,417,452,484, 681]]))
major_mutations = mutations.apply(lambda x:','.join([y for y in x.split(',') if y and int(y[1:-1]) in [154,681]]))
meta_muts = pd.concat([major_mutations, meta.loc[:, ["date", "date_submitted", "country", "region", "Nextstrain_clade", "pango_lineage"]]], axis=1).loc[mutations.index]

for muts, group in meta_muts.groupby(by='S'):
    if len(group)>5:
        print(f"{muts:<20}\t{len(group)}\t{[(x,c) for x,c in sorted(Counter(group['pango_lineage']).items(), key=lambda x:-x[1])[:10] if c>3]}")
