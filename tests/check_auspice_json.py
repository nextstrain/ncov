import argparse
import json
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Ensure certain values are present for a given node trait",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--json', type=str, metavar="JSON", required=True, help="Auspice JSON")
    parser.add_argument('--attr', type=str, metavar="KEY", required=True, help="node attr to collect")
    parser.add_argument('--values', type=str, nargs="+", metavar="VALUE", required=True, help="values to check")
    args = parser.parse_args()

    values_seen = set()
    
    def collect(node):
        v = node.get("node_attrs", {}).get(args.attr, {}).get("value", "")
        if v:
            values_seen.add(v)
        for child in node.get("children", []):
            collect(child)

    with open(args.json, "r") as f:
        input_json = json.load(f)

    collect(input_json["tree"])

    if not values_seen >= set(args.values):
        print("Following values missing from JSON:", set(args.values)-values_seen)
        sys.exit(1)
