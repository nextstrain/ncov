#!/usr/bin/env python3
import argparse
from augur.utils import read_tree, read_node_data, read_metadata, write_json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find polytomies in a given tree that all belong to the same metadata group",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--metadata", required=True, help="metadata")
    parser.add_argument("--mutations", required=True, help="mutations node data JSON")
    parser.add_argument("--attribute-name", default="cluster_id", help="name of attribute to store in output JSON")
    parser.add_argument("--group-by", default="division", help="identify polytomies where all tips are in the same group")
    parser.add_argument("--min-tips", type=int, default=3, help="minimum tips per polytomy to be consider as a cluster")
    parser.add_argument("--output-node-data", required=True, help="node data JSON with cluster id")

    args = parser.parse_args()

    tree = read_tree(args.tree)
    tree.collapse_all(lambda c: c.branch_length < 1e-5)

    metadata, columns = read_metadata(args.metadata)
    muts = read_node_data(args.mutations)
    attribute_name = args.attribute_name
    group_by = args.group_by

    polytomies = []
    for node in tree.find_clades(terminal=False):
        if node == tree.root:
            continue

        any_muts = False
        groups = set()
        children = 0
        for child in node.clades:
            if child.is_terminal() and child.name:
                children += 1
                any_muts |= (len(muts["nodes"].get(child.name, {}).get("muts", [])) > 0)
                groups.add(metadata[child.name][group_by])

        if not any_muts and children >= args.min_tips and len(groups) == 1:
            polytomies.append(node)

    node_data = {}
    clusters = 0
    for polytomy in polytomies:
        if polytomy.name:
            node_data[polytomy.name] = {
                attribute_name: f"cluster_{clusters}",
            }

        for child in polytomy.clades:
            if child.is_terminal():
                node_data[child.name] = {
                    attribute_name: f"cluster_{clusters}",
                }

        clusters += 1

    write_json({"nodes": node_data}, args.output_node_data)
