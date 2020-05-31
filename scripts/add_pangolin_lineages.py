import argparse
import json
import requests
import io
import pandas as pd
from Bio import Phylo

URL = "https://raw.githubusercontent.com/hCoV-2019/lineages/master/lineages/data/lineages.metadata.csv"
def get_pangolin_lineages():
    r  = requests.get(URL)
    if not r.ok:
        print(f"Failed to fetch {URL}", file=sys.stderr)
        exit(1)
        r.close()

    d = pd.read_csv(io.StringIO(r.text))
    lineages = {}
    for ri, row in d.iterrows():
        if row['lineage']:
            lineages[row['name']] = row['lineage']
            lineages[row['name'].replace('_', '')] = row['lineage']
    return lineages

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add pangolin lineage definition to tree json as coloring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="tree file")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    lineages =get_pangolin_lineages()

    T = Phylo.read(args.tree, 'newick')
    node_data = {}
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            if n.name in lineages:
                n.pangolin = lineages[n.name]
            else:
                n.pangolin = 'unassigned'
        else:
            daughter_clades = list({c.pangolin for c in n if c.pangolin!='unassigned'})
            if len(daughter_clades)==1:
                n.pangolin=daughter_clades[0]
            else:
                n.pangolin = 'ambiguous'
        if n.pangolin!='ambiguous':
            node_data[n.name] = {'pangolin_lineage':n.pangolin}


    with open(args.output, 'w', encoding='utf-8') as f:
        json.dump({"nodes":node_data}, f, indent=2)
