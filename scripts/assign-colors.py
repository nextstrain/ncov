import argparse
import pandas as pd

# Forced colours MUST NOT appear in the ordering TSV
forced_colors = {
  "division": {
    "Diamond Princess": "#CCCCCC",
  },
  "location": {
    "Diamond Princess": "#CCCCCC",
  },
  "division": {
    "Grand Princess": "#AAAAAA",
  },
  "location": {
    "Grand Princess": "#AAAAAA",
  }
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="input ordering file")
    parser.add_argument('--color-schemes', type=str, required=True, help="input color schemes file")
    parser.add_argument('--output', type=str, required=True, help="output colors tsv")
    parser.add_argument('--metadata', type=str, help="if provided, restrict colors to only those found in metadata")
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

    if(args.metadata):
        metadata=pd.read_csv(args.metadata, delimiter='\t')
        for k in assignment.keys():
            if k in metadata:
                new_list = [x for x in assignment[k] if x in metadata[k].unique()]
                assignment[k] = new_list

        if 'focal_country' in metadata:
            assignment['focal_country'] = assignment['country']

    schemes = {}
    counter = 0
    with open(args.color_schemes) as f:
        for line in f.readlines():
            counter += 1
            array = line.lstrip().rstrip().split("\t")
            schemes[counter] = array

    with open(args.output, 'w') as f:
        for trait_name, trait_array in assignment.items():
            color_array = schemes[len(trait_array)]
            extra_trait_values = list(forced_colors.get(trait_name, {}).keys())
            extra_color_values = list(forced_colors.get(trait_name, {}).values())

            zipped = list(zip(trait_array+extra_trait_values, color_array+extra_color_values))
            for trait_value, color in zipped:
                f.write(trait_name + "\t" + trait_value + "\t" + color + "\n")
            f.write("\n")
