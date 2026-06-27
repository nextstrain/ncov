import argparse
import json
import pandas as pd

import color_utils

# Forced colours MUST NOT appear in the ordering TSV
forced_colors = {
}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="input ordering file")
    parser.add_argument('--color-schemes', type=str, required=True, help="input color schemes file")
    parser.add_argument('--metadata', type=str, help="if provided, restrict colors to only those found in metadata")
    parser.add_argument('--clade-node-data', type=str, help="if provided, restrict to only those clades found in tree")
    parser.add_argument('--clade-recency', type=color_utils.relative_date, metavar='DURATION',
        help="""if provided, restrict to clades found in tree within this time
             frame. Format: ISO 8601 duration with optional P prefix (e.g. '1W',
             'P1W')""")
    parser.add_argument('--output', type=str, required=True, help="output colors tsv")
    args = parser.parse_args()

    assignment = {}
    with open(args.ordering) as f:
        for line in f.readlines():
            array = line.lstrip().rstrip().split("\t")
            if len(array) == 2:
                name = array[0]
                trait = array[1]
                if name not in assignment:
                    assignment[name] = [trait]
                else:
                    assignment[name].append(trait)

    # if metadata supplied, go through and
    # 1. remove assignments that don't exist in metadata
    # 2. remove assignments that have 'focal' set to 'False' in metadata
    if args.metadata:
        metadata = pd.read_csv(args.metadata, delimiter='\t')
        for name, trait in assignment.items():
            if name in metadata['strain'].values:
                subset_present = [x for x in assignment[name] if x in metadata[name].unique()]
                assignment[name] = subset_present
            if name in metadata['strain'].values and 'focal' in metadata.columns:
                focal_list = metadata.loc[metadata['focal'] == True, name].unique()
                subset_focal = [x for x in assignment[name] if x in focal_list]
                assignment[name] = subset_focal

    # if node json is supplied, restrict clade_membership to clades in the tree within
    # the specified recency (shared with assign-lineage-colors.py via color_utils)
    if args.clade_node_data and "clade_membership" in assignment:
        with open(args.clade_node_data) as fh:
            clades = json.load(fh)['nodes']
        metadata = pd.read_csv(args.metadata, delimiter='\t') if args.metadata else None
        subset_present = color_utils.present_clades(clades, metadata, args.clade_recency)
        assignment["clade_membership"] = [x for x in assignment["clade_membership"]
                                          if x in subset_present]

    schemes = color_utils.load_color_schemes(args.color_schemes)

    with open(args.output, 'w') as f:
        for trait_name, trait_array in assignment.items():
            if len(trait_array) == 0:
                print(f"No traits found for {trait_name}")
                continue
            if len(schemes) < len(trait_array):
                print(f"WARNING: insufficient colours available for trait {trait_name} - reusing colours!")
            color_array = color_utils.assign_scheme(trait_array, schemes)
            extra_trait_values = list(forced_colors.get(trait_name, {}).keys())
            extra_color_values = list(forced_colors.get(trait_name, {}).values())

            zipped = list(zip(trait_array+extra_trait_values, color_array+extra_color_values))
            for trait_value, color in zipped:
                f.write(trait_name + "\t" + trait_value + "\t" + color + "\n")
            f.write("\n")
