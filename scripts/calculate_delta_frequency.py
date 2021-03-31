"""Calculate the change in frequency for clades over time (aka the delta frequency or dfreq).
"""
import argparse
from augur.frequency_estimators import TreeKdeFrequencies
from augur.utils import read_node_data, read_tree, write_json
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np


def read_frequencies(frequencies_file):
    """Returns a dictionary of frequencies and their parameters indexed by strain
    name from a given auspice tip frequencies file.

    """
    with open(frequencies_file) as fh:
        frequencies_json = json.load(fh)

    parameters = {}
    frequencies = {}

    for key, values in frequencies_json.items():
        if "frequencies" in values:
            frequencies[key] = values["frequencies"]
        else:
            parameters[key] = values

    return frequencies, parameters


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the change in frequency for clades over time",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick tree")
    parser.add_argument("--frequencies", required=True, help="frequencies JSON")
    parser.add_argument("--delta-pivots", type=int, default=1, help="number of frequency pivots to look back in time for change in frequency calculation")
    parser.add_argument("--attribute-name", default="delta_frequency", help="name of the annotation to store in the node data JSON output")
    parser.add_argument("--include-tips", action="store_true", help="include change of frequency for tips in output. This output tends to be less meaningful than change of frequency for internal nodes (i.e., clades).")
    parser.add_argument("--output", required=True, help="JSON of delta frequency annotations for nodes in the given tree")

    args = parser.parse_args()

    # Load the tree.
    tree = read_tree(args.tree)

    # Load frequencies.
    frequencies, parameters = read_frequencies(args.frequencies)
    pivots = parameters["pivots"]

    # Determine the total time that elapsed between the current and past timepoint.
    delta_time = pivots[-1] - pivots[-(args.delta_pivots + 1)]

    # Calculate frequencies for internal nodes by summing the frequencies of
    # their respective tips. Then calculate the change in frequency for each
    # node from the resulting frequencies.
    delta_frequency = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            # We already know the frequencies of each terminal node, so
            # store those frequencies with the corresponding node of the tree.
            node.frequencies = frequencies[node.name]
        else:
            # For each internal node, sum the frequencies of its immediate
            # children. Since we are walking through the tree from the bottom
            # up, each child node will always have frequencies calculated
            # before its parent node. Thus, summing the frequencies of the
            # immediate children in postorder gives us the sum of the frequencies
            # of all children of a node (not just the immediate ones).
            node.frequencies = np.array([
                child.frequencies
                for child in node.clades
            ]).sum(axis=0)

        if not node.is_terminal() or args.include_tips:
            # Calculate the change in frequency over the requested time period.
            node_delta_frequency = (node.frequencies[-1] - node.frequencies[-(args.delta_pivots + 1)]) / delta_time

            if node_delta_frequency != 0:
                delta_frequency[node.name] = {
                    args.attribute_name: node_delta_frequency,
                    "current_frequency": node.frequencies[-1]
                }

    # Write out the node annotations.
    write_json({"nodes": delta_frequency}, args.output)
