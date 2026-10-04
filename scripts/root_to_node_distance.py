#!/usr/bin/env python3
"""
Sum per-branch sequence distances along the path from the root to each node.

`augur distance --compare-to root` compares every node directly against the root
sequence, so recurrent and reverting substitutions at the same site do not
accumulate -- the metric saturates as sites are hit repeatedly. This script
instead walks the tree and, for every edge, adds the per-branch distance
computed with augur's own distance-map logic, so each substitution event along
the root-to-node path is counted. The result is a molecular-clock-like tally of
cumulative amino-acid changes rather than a net difference from the reference.
"""
import argparse
from collections import defaultdict

from Bio import Phylo
from augur.distance import read_distance_map, get_distance_between_nodes
from augur.reconstruct_sequences import load_alignments
from augur.io import write_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sum per-branch sequence distances from the root to each node.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--tree", required=True, help="Newick tree with named internal nodes")
    parser.add_argument("--alignment", nargs="+", required=True,
                        help="per-gene FASTA alignment(s) including internal-node sequences")
    parser.add_argument("--gene-names", nargs="+", required=True,
                        help="gene name for each alignment, paired positionally with --alignment")
    parser.add_argument("--map", required=True,
                        help="distance map JSON (sites, weights, ignored_characters)")
    parser.add_argument("--attribute-name", required=True,
                        help="name to store the cumulative distance under in the node-data JSON")
    parser.add_argument("--output", required=True, help="output node-data JSON")
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    distance_map = read_distance_map(args.map)

    # Flatten per-gene alignments to {node_name: {gene: sequence}}, exactly as augur distance does.
    alignments = load_alignments(args.alignment, args.gene_names)
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Walk from the root, accumulating the per-branch distance on each edge.
    # Iterative preorder guarantees a parent is assigned before its children and,
    # unlike recursion, is safe for the deep trees ncov can produce.
    cumulative = {tree.root.name: 0}
    node_data = {tree.root.name: {args.attribute_name: 0}}
    for node in tree.find_clades(order="preorder"):
        for child in node.clades:
            branch_distance = get_distance_between_nodes(
                sequences_by_node_and_gene[node.name],
                sequences_by_node_and_gene[child.name],
                distance_map,
            )
            cumulative[child.name] = cumulative[node.name] + branch_distance
            node_data[child.name] = {args.attribute_name: cumulative[child.name]}

    write_json({"nodes": node_data}, args.output)
