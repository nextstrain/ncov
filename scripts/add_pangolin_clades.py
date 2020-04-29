import argparse
import json


def add_clade(c, clades):
    if c["name"] in clades:
        c["node_attrs"]["pangolin-clade"] = {"value": clades[c["name"]]}
    if "children" in c:
        for n in c["children"]:
            add_clade(n, clades)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--clades', type=str, required=True, help="clade file")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    clades = {}
    with open(args.clades, 'r') as fh:
        for line in fh:
            strain, clade, _1, _2, status  = line.strip().split(",")[:5]
            if status=='success':
                clades[strain]=clade

    with open(args.input, "r") as f:
        input_json = json.load(f)

    keys_to_remove = ["genbank_accession", "gisaid_epi_isl"]

    fixed_colorings = []
    for coloring in input_json["meta"]["colorings"]:
        if coloring['key'] not in keys_to_remove:
            fixed_colorings.append(coloring)

    input_json["meta"]["colorings"] = fixed_colorings

    add_clade(input_json["tree"], clades)
    input_json["meta"]["colorings"].append(
      {
        "key": "pangolin-clade",
        "title": "pangolin-clade",
        "type": "categorical"
      })

    with open(args.output, 'w') as f:
        json.dump(input_json, f, indent=2)
