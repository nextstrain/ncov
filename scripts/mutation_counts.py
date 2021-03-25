import argparse, json
from Bio import Phylo, SeqIO
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="count distances based on distance maps with counting gaps as one event",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="tree file")
    parser.add_argument('--alignment', type=str)
    parser.add_argument('--gene-names', nargs='+')
    parser.add_argument('--compare-to', nargs='+')
    parser.add_argument('--map', nargs='+')
    parser.add_argument('--attribute-name', nargs='+')
    parser.add_argument('--output', type=str, required=True, help="output json file")
    args = parser.parse_args()

    if type(args.map)==str:
        args.map = [args.map]

    if type(args.attribute_name)==str:
        args.attribute_name = [args.attribute_name]
    
    node_data = {}

    T = Phylo.read(args.tree, 'newick')
    aln = {s.name:str(s.seq) for s in SeqIO.parse(args.alignment, 'fasta')}
    root = aln[T.root.name]

    for map_file, attribute in zip(args.map, args.attribute_name):
        with open(map_file) as fh:
            map = json.load(fh)

        weights = {int(k)-1:v for k,v in map["map"][args.gene_names[0]].items()}
        default = map['default']
        for n in T.find_clades():
            if n.name not in node_data: node_data[n.name] = {}
            distance = 0
            gap_extend=False
            if n.name not in aln:
                print(f"{n.name} not found in alignment")
                continue

            for p, (a,d) in enumerate(zip(root, aln[n.name])):
                if a!=d and (not gap_extend):
                    distance+=weights.get(p, default)

                gap_extend= d=='-'
            node_data[n.name][f"{attribute}"] = distance

    with open(args.output, 'w') as fh:
        json.dump({"nodes":node_data}, fh)
