import argparse, json
from Bio import Phylo
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="transform nextalign output to sparse format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="tree file")
    parser.add_argument('--mutation-summary')
    parser.add_argument('--genes', nargs='+')
    parser.add_argument('--output', type=str, required=True, help="output json file")
    args = parser.parse_args()

    mutations = pd.read_csv(args.mutation_summary, sep='\t', index_col=0).fillna('')

    d = {}
    T = Phylo.read(args.tree, 'newick')
    for n in T.get_terminals():
        d[n.name] = {}
        try:
            m = mutations.loc[n.name]
        except:
            print(f"no data found for {n.name}")
            continue
        for gene in args.genes:
            d[n.name][f"{gene}_mutation_count"] = len([x for x in m[gene].split(',') if x and x[-1]!='-'])

    with open(args.output, 'w') as fh:
        json.dump({"nodes":d}, fh)
