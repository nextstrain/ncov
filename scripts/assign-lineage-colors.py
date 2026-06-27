"""
Assign colors to Pango lineages (Nextclade_pango, pango_lineage) so they correspond to,
and are ordered by, the Nextstrain clade coloring.

A lineage is colored if and only if its clade is colored by clade_membership — i.e. the
clade is in the clade_membership ordering AND has a sequence within --clade-recency (the
exact same set scripts/assign-clade-colors.py colors). Lineages are grouped by clade (in the
clade_membership order) and, within a clade, ordered by the Nextclade reference tree's
pre-order (from the build's Nextclade dataset zip), which groups sublineages — including
Pango aliases (e.g. XFG.3.4.1 and its alias QF.2) — adjacently, so colors follow clade
structure. Output rows are appended to the build's colors.tsv.
"""

import argparse
import json
import zipfile
from collections import Counter, defaultdict

import pandas as pd

import color_utils

TRAITS = ["Nextclade_pango", "pango_lineage"]


def clade_membership_order(ordering_path):
    """The clade_membership values from color_ordering.tsv, in order (the colorable clades)."""
    order = []
    with open(ordering_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0] == "clade_membership":
                order.append(parts[1])
    return order


def lineage_preorder(nextclade_dataset_zip):
    """lineage -> pre-order (DFS) index from the dataset's reference tree.json."""
    tree = json.loads(zipfile.ZipFile(nextclade_dataset_zip).read("tree.json"))["tree"]
    order, counter = {}, [0]

    def walk(node):
        attr = node.get("node_attrs", {}).get("Nextclade_pango")
        lineage = attr.get("value") if isinstance(attr, dict) else None
        if lineage and lineage not in order:
            order[lineage] = counter[0]
            counter[0] += 1
        for child in node.get("children", []):
            walk(child)

    walk(tree)
    return order


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clade-node-data", required=True, help="clades.json node data (strain -> clade_membership)")
    p.add_argument("--metadata", required=True, help="build metadata (strain, date, Nextclade_pango, pango_lineage)")
    p.add_argument("--clade-recency", type=color_utils.relative_date, metavar="DURATION",
                   help="ISO 8601 duration; restrict to clades with a sequence in this window")
    p.add_argument("--color-schemes", required=True, help="color schemes TSV")
    p.add_argument("--nextclade-dataset", required=True, help="Nextclade dataset zip (for tree.json structure)")
    p.add_argument("--clade-ordering", required=True, help="color_ordering.tsv (for clade_membership order)")
    p.add_argument("--input-colors", help="existing colors TSV (clade/geographic) to prepend to the output")
    p.add_argument("--output", required=True, help="output colors TSV (input-colors plus Pango lineage colors)")
    args = p.parse_args()

    with open(args.clade_node_data) as fh:
        clade_nodes = json.load(fh)["nodes"]
    strain_clade = {s: info["clade_membership"] for s, info in clade_nodes.items()
                    if info.get("clade_membership")}
    metadata = pd.read_csv(args.metadata, delimiter="\t")

    # Colored clades = the clade_membership ordering restricted to clades present within
    # the recency window — exactly the set assign-clade-colors.py colors.
    present = color_utils.present_clades(clade_nodes, metadata, args.clade_recency)
    clade_order = [c for c in clade_membership_order(args.clade_ordering) if c in present]
    clade_rank = {c: i for i, c in enumerate(clade_order)}

    preorder = lineage_preorder(args.nextclade_dataset)
    schemes = color_utils.load_color_schemes(args.color_schemes)

    with open(args.output, "w") as out:
        if args.input_colors:
            with open(args.input_colors) as base:
                out.write(base.read())
        for trait in TRAITS:
            if trait not in metadata.columns:
                continue
            # present lineage -> its clade (majority clade across the lineage's strains)
            votes = defaultdict(Counter)
            for strain, lineage in zip(metadata["strain"], metadata[trait]):
                clade = strain_clade.get(strain)
                if clade and pd.notna(lineage):
                    votes[lineage][clade] += 1
            lineage_clade = {lin: c.most_common(1)[0][0] for lin, c in votes.items()}

            colored = [lin for lin, clade in lineage_clade.items() if clade in clade_rank]
            # group by clade (clade order), then ref-tree pre-order, then name
            colored.sort(key=lambda lin: (clade_rank[lineage_clade[lin]],
                                          preorder.get(lin, 10 ** 9), lin))

            colors = color_utils.assign_scheme(colored, schemes)
            for lineage, color in zip(colored, colors):
                out.write(f"{trait}\t{lineage}\t{color}\n")
            out.write("\n")


if __name__ == "__main__":
    main()
