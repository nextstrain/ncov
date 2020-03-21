import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        input_json = json.load(f)

    keys_to_remove = ["genbank_accession", "gisaid_epi_isl"]

    fixed_colorings = []
    for coloring in input_json["meta"]["colorings"]:
        if coloring['key'] not in keys_to_remove:
            fixed_colorings.append(coloring)

    input_json["meta"]["colorings"] = fixed_colorings

    with open(args.output, 'w') as f:
        json.dump(input_json, f, indent=2)
