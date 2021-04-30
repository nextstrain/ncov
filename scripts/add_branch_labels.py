import argparse
import json
from Bio import Phylo
from collections import defaultdict

def extract_spike_mutations(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        smuts = node.get("aa_muts", {}).get("S", [])
        if smuts:
            data[name] = ", ".join(smuts)
    return data

def extract_clade_labels(node_data):
    data = {}
    for name, node in node_data["nodes"].items():
        if "clade_annotation" in node:
            data[name] = node["clade_annotation"]
    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--mutations', type=str, required=False, help="mutations node data file")
    parser.add_argument('--emerging-clades', type=str, required=True, help="emerging clades node data file")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        auspice_json = json.load(f)

    if args.mutations:
        with open(args.mutations, "r") as f:
            spike_mutations = extract_spike_mutations(json.load(f))
    else:
        spike_mutations = {}

    with open(args.emerging_clades, "r") as f:
        clade_labels = extract_clade_labels(json.load(f))

    def attach_labels(n): # closure
      if n["name"] in spike_mutations or n["name"] in clade_labels:
          if "branch_attrs" not in n:
              n["branch_attrs"]={}
          if "labels" not in n["branch_attrs"]:
              n["branch_attrs"]["labels"]={}
          if n["name"] in spike_mutations:
              n["branch_attrs"]["labels"]["spike_mutations"] = spike_mutations[n["name"]]
          if n["name"] in clade_labels:
              n["branch_attrs"]["labels"]["emerging_lineage"] = clade_labels[n["name"]]

      if "children" in n:
          for c in n["children"]:
              attach_labels(c)

    attach_labels(auspice_json["tree"])

    with open(args.output, 'w') as f:
        json.dump(auspice_json, f, indent=2)
