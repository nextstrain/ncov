"""Calculate the change in frequency for clades over time (aka the delta frequency or dfreq).
"""
import argparse
from augur.frequency_estimators import logit_transform
from augur.utils import annotate_parents_for_tree, read_node_data, read_tree, write_json
import Bio.Phylo
from collections import defaultdict
import json
import numpy as np
from scipy.stats import linregress
import sys


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
    parser.add_argument(
        "--method",
        default="linear",
        choices=("linear", "logistic"),
        help="""method to use when calculating slope of frequency changes per clade.
        The 'linear' method calculates the slope between the most recent timepoint and the timepoint associted with the number of pivots back in time requested by `--delta-pivots`.
        The 'logistic' method fits applies logistic regression per clade to frequencies of each timepoint between the latest and earliest requested timepoint and uses the slope from this regression."""
    )
    parser.add_argument(
        "--min-tips",
        default=10,
        type=int,
        help="minimum number of tips for internal nodes to calculate delta frequency for. Nodes below this number inherit the values of their parent node."
    )
    parser.add_argument("--attribute-name", default="delta_frequency", help="name of the annotation to store in the node data JSON output")
    parser.add_argument("--include-tips", action="store_true", help="include change of frequency for tips in output. This output tends to be less meaningful than change of frequency for internal nodes (i.e., clades).")
    parser.add_argument("--output", required=True, help="JSON of delta frequency annotations for nodes in the given tree")

    args = parser.parse_args()

    # Load the tree.
    tree = read_tree(args.tree)
    tree = annotate_parents_for_tree(tree)

    # Load frequencies.
    frequencies, parameters = read_frequencies(args.frequencies)
    pivots = np.array(parameters["pivots"])

    # Determine the total time that elapsed between the current and past timepoint.
    first_pivot_index = -(args.delta_pivots + 1)
    last_pivot_index = -1
    delta_time = pivots[last_pivot_index] - pivots[first_pivot_index]

    # Calculate frequencies for internal nodes by summing the frequencies of
    # their respective tips.
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            # We already know the frequencies of each terminal node, so
            # store those frequencies with the corresponding node of the tree.
            node.frequencies = np.array(frequencies[node.name])
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

    # Calculate the change in frequency for each node from the precalculated
    # frequencies. The signal for smaller clades is noisier, so we set a minimum
    # clade frequency below which clades/tips inherit their parent's delta
    # frequency value.
    delta_frequency = {}
    for node in tree.find_clades(order="preorder"):
        # Always annotate the current frequency of each node.
        delta_frequency[node.name] = {
            "current_frequency": node.frequencies[last_pivot_index]
        }

        if node.count_terminals() >= args.min_tips:
            # Calculate the change in frequency over the requested time period.
            if args.method == "linear":
                node_delta_frequency = (node.frequencies[last_pivot_index] - node.frequencies[first_pivot_index]) / delta_time
            elif args.method == "logistic":
                x_pivots = pivots[first_pivot_index:]

                # Transform most recent frequencies prior to fitting linear
                # regression to better represent logistic growth we expect from
                # SARS-CoV-2 clades. This transformation accounts for numerical
                # error with its second argument to avoid infinite values in the
                # transform (as when frequencies equal 0 or 1).
                y_frequencies = logit_transform(
                    node.frequencies[first_pivot_index:],
                    pc=1e-4
                )

                # Fit linear regression to pivots and frequencies and use the
                # resulting slope as the measure of recent clade growth or
                # decline.
                model = linregress(x_pivots, y_frequencies)
                node_delta_frequency = model.slope
            else:
                print(f"Error: The request method, '{args.method}', is not supported.", file=sys.stderr)
                sys.exit(1)

            delta_frequency[node.name][args.attribute_name] = node_delta_frequency
        else:
            # If the current node is low frequency, use its parent node's delta frequency value.
            delta_frequency[node.name][args.attribute_name] = delta_frequency[node.parent.name][args.attribute_name]

    # Write out the node annotations.
    write_json({"nodes": delta_frequency}, args.output)