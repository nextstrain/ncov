import argparse
import datetime
import isodate
import pandas as pd

# Forced colours MUST NOT appear in the ordering TSV
forced_colors = {
}

def date_within_last_n_months(date_str, cutoff_date):
    if 'XX' in date_str:
        return False  # Ignore uncertain dates
    try:
        date = datetime.datetime.strptime(date_str, "%Y-%m-%d").date()
        return date >= cutoff_date
    except ValueError:
        return False


def relative_date(duration: str):
    """
    Convert an ISO 8601 duration to an absolute date by subtracting it from the
    current date.

    `duration` should be a backwards-looking relative date in ISO 8601 duration
    format with optional P prefix (e.g. '1W', 'P1W').
    """
    if duration.startswith('P'):
        duration = duration
    else:
        duration = 'P' + duration
    return datetime.date.today() - isodate.parse_duration(duration)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="input ordering file")
    parser.add_argument('--color-schemes', type=str, required=True, help="input color schemes file")
    parser.add_argument('--metadata', type=str, help="if provided, restrict colors to only those found in metadata")
    parser.add_argument('--clade-node-data', type=str, help="if provided, restrict to only those clades found in tree")
    parser.add_argument('--clade-recency', type=relative_date, metavar='DURATION',
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

    # if node json is supplied, restrict to clades names in the tree within the specified recency
    if args.clade_node_data and "clade_membership" in assignment:
        with open(args.clade_node_data) as fh:
            import json
            clades = json.load(fh)['nodes']

        if args.clade_recency is not None and args.metadata:
            # Generate a set of present values within the specified recency
            subset_present = set()
            metadata = pd.read_csv(args.metadata, delimiter='\t')
            for strain, info in clades.items():
                if strain in metadata['strain'].values:
                    date_str = metadata.loc[metadata['strain'] == strain, 'date'].values[0]
                    if date_within_last_n_months(date_str, args.clade_recency):
                        subset_present.add(info["clade_membership"])

            # Restrict to only those present while maintaining order
            assignment["clade_membership"] = [x for x in assignment["clade_membership"]
                                              if x in subset_present]
        else:
            # If no clade_recency is provided, look for all clades present in the tree
            subset_present = set([x["clade_membership"] for x in clades.values()])
            assignment["clade_membership"] = [x for x in assignment["clade_membership"]
                                              if x in subset_present]

    schemes = {}
    counter = 0
    with open(args.color_schemes) as f:
        for line in f.readlines():
            counter += 1
            array = line.lstrip().rstrip().split("\t")
            schemes[counter] = array

    with open(args.output, 'w') as f:
        for trait_name, trait_array in assignment.items():
            if len(trait_array)==0:
                print(f"No traits found for {trait_name}")
                continue
            if len(schemes)<len(trait_array):
              print(f"WARNING: insufficient colours available for trait {trait_name} - reusing colours!")
              remain = len(trait_array)
              color_array = []
              while(remain>0):
                if (remain>len(schemes)):
                  color_array = [*color_array, *schemes[len(schemes)]]
                  remain -= len(schemes)
                else:
                  color_array = [*color_array, *schemes[remain]]
                  remain = 0
            else:
              color_array = schemes[len(trait_array)]
            extra_trait_values = list(forced_colors.get(trait_name, {}).keys())
            extra_color_values = list(forced_colors.get(trait_name, {}).values())

            zipped = list(zip(trait_array+extra_trait_values, color_array+extra_color_values))
            for trait_value, color in zipped:
                f.write(trait_name + "\t" + trait_value + "\t" + color + "\n")
            f.write("\n")
