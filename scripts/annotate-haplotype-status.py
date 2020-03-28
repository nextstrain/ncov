"""Annotate haplotype status for all sequences that match a given reference node name.
"""
import argparse
from augur.utils import write_json
import json
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ancestral-sequences", required=True, help="node data JSON of nucleotide mutations including observed and inferred by TreeTime")
    parser.add_argument("--reference-node-name", required=True, help="name of the node whose sequence flags the desired haplotype")
    parser.add_argument("--attribute-name", default="haplotype_status", help="name of attribute for haplotype status")
    parser.add_argument("--output", required=True, help="node data JSON with annotated haplotype status based on the given reference node's sequence")

    args = parser.parse_args()

    with open(args.ancestral_sequences, "r") as fh:
        sequences = json.load(fh)

    if args.reference_node_name not in sequences["nodes"]:
        print("ERROR: Could not find the requested reference node named '%s' in the given ancestral sequences." % args.reference_node_name, file=sys.stderr)
        sys.exit(1)

    haplotype_sequence = sequences["nodes"][args.reference_node_name]["sequence"]
    haplotype_status = {"nodes": {}}

    for node in sequences["nodes"]:
        if sequences["nodes"][node]["sequence"] == haplotype_sequence:
            status = "haplotype matches %s" % args.reference_node_name
        else:
            status = "haplotype does not match %s" % args.reference_node_name

        haplotype_status["nodes"][node] = {args.attribute_name: status}

    write_json(haplotype_status, args.output)
