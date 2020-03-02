from Bio import Phylo
import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate TMRCA assuming a star topology and a poisson mutation process",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick file (only topology is used)")
    parser.add_argument("--node-data", required=True, help="JSON from `augur refine`")
    parser.add_argument("--output", required=True, help="Newick tree file output")
    args = parser.parse_args()


    with open(args.node_data) as fh:
        j = json.load(fh)
    
    T = Phylo.read(args.tree, 'newick')

    for n in T.find_clades():
        n.branch_length = j["nodes"][n.name]["mutation_length"]
        print(n.name, n.branch_length)

    Phylo.write(T, args.output, "newick")