#!/usr/bin/env python3
import argparse
from augur.io import read_metadata
from augur.utils import read_tree, read_node_data
from collections import Counter
import csv
import hashlib

MAX_HASH_LENGTH = 7


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
    parser.add_argument("--output", required=True, help="tab-delimited file with strain, cluster id, and group value for each strain")

    args = parser.parse_args()

    tree = read_tree(args.tree)
    tree.collapse_all(lambda c: c.branch_length < 1e-5)

    metadata = read_metadata(args.metadata)
    muts = read_node_data(args.mutations)
    attribute_name = args.attribute_name
    group_by = args.group_by

    polytomies = []
    for node in tree.find_clades(terminal=False):
        if node == tree.root:
            continue

        count_by_group = Counter()
        polytomy_sequence_id = None
        for child in node.clades:
            if child.is_terminal() and child.name:
                child_muts_data = muts["nodes"].get(child.name, {})
                any_muts = (len(child_muts_data.get("muts", [])) > 0)
                if not any_muts:
                    count_by_group[metadata.loc[child.name, group_by]] += 1

                    if polytomy_sequence_id is None and "sequence" in child_muts_data:
                        polytomy_sequence_id = hashlib.sha256(child_muts_data["sequence"].encode()).hexdigest()[:MAX_HASH_LENGTH]

        if any(count >= args.min_tips for count in count_by_group.values()):
            polytomies.append({"node": node, "name": polytomy_sequence_id})

    with open(args.output, "w") as oh:
        writer = csv.DictWriter(
            oh,
            fieldnames=(
                "strain",
                args.attribute_name,
                group_by
            ),
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        clusters = 0
        for polytomy_data in polytomies:
            polytomy = polytomy_data["node"]
            polytomy_sequence_id = polytomy_data["name"]

            if polytomy.name:
                writer.writerow({
                    "strain": polytomy.name,
                    args.attribute_name: polytomy_sequence_id,
                    group_by: metadata.loc[polytomy.name, group_by]
                })

            for child in polytomy.clades:
                if child.is_terminal():
                    writer.writerow({
                        "strain": child.name,
                        args.attribute_name: polytomy_sequence_id,
                        group_by: metadata.loc[child.name, group_by]
                    })

            clusters += 1
