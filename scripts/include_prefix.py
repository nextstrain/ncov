import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rename strains to include specified prefix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-auspice', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--input-tip-frequencies', type=str, metavar="JSON", required=True, help="input tip frequencies JSON")
    parser.add_argument("--prefix", type=str, nargs='?', const='', help="prefix to add to strains")
    parser.add_argument('--output-auspice', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    parser.add_argument('--output-tip-frequencies', type=str, metavar="JSON", required=True, help="output tip frequencies JSON")
    args = parser.parse_args()

    # update Auspice JSON
    with open(args.input_auspice, "r") as f:
        auspice_json = json.load(f)

    if args.prefix:
        def update_strain_names(n): # closure
            if "NODE_" not in n["name"] and args.prefix not in n["name"]:
                n["name"] = args.prefix + n["name"]

            if "children" in n:
                for c in n["children"]:
                    update_strain_names(c)
        update_strain_names(auspice_json["tree"])

    with open(args.output_auspice, 'w') as f:
        json.dump(auspice_json, f, indent=2)

    # update tip frequencies JSON
    with open(args.input_tip_frequencies, "r") as f:
        tip_frequencies_json = json.load(f)

    if args.prefix:
        modified_tip_frequencies_json = {}
        for key in tip_frequencies_json:
            if key != "generated_by" and key != "pivots":
                if "NODE_" not in key and args.prefix not in key:
                    modified_tip_frequencies_json[args.prefix + key] = tip_frequencies_json[key]
                else:
                    modified_tip_frequencies_json[key] = tip_frequencies_json[key]
            else:
                modified_tip_frequencies_json[key] = tip_frequencies_json[key]
    else:
        modified_tip_frequencies_json = tip_frequencies_json

    with open(args.output_tip_frequencies, 'w') as f:
        json.dump(modified_tip_frequencies_json, f, indent=2)
