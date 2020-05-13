import argparse
import json


def add_lineage(c, clades):
    if c["name"] in clades:
        c["node_attrs"]["pangolin-lineage"] = {"value": clades[c["name"]]}
    if "children" in c:
        for n in c["children"]:
            add_lineage(n, clades)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add pangolin lineage definition to tree json as coloring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--lineages', type=str, required=True, help="csv file with pangolin lineages")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    # read file with pangolin lineage definitions
    lineages = {}
    with open(args.lineages, 'r') as fh:
        for line in fh:
            strain, country, travel_history, sample_date, epiweek, lineage, representative = line.strip().split(",")
            lineages[strain] = lineage

    with open(args.input, "r") as f:
        input_json = json.load(f)

    add_lineage(input_json["tree"], lineages)
    input_json["meta"]["colorings"].append(
      {
        "key": "pangolin-lineage",
        "title": "pangolin-lineage",
        "type": "categorical"
      })

    with open(args.output, 'w') as f:
        json.dump(input_json, f, indent=2)
