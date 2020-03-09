from collections import defaultdict
from itertools import chain
import argparse

# Forced colors MUST NOT appear in the ordering TSV
forced_colors = {
  "division": {
    "Diamond Princess": "#CCCCCC",
  },
  "location": {
    "Diamond Princess": "#CCCCCC",
  }
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', required=True,
                        help="input ordering file")
    parser.add_argument('--color-schemes', required=True,
                        help="input color schemes file")
    parser.add_argument('--output', required=True,
                        help="output colors tsv")
    args = parser.parse_args()

    assignment = defaultdict(list)
    with open(args.ordering) as f:
        for line in f:
            try:
                name, trait = line.strip().split("\t")
            except ValueError:
                pass
            else:
                assignment[name].append(trait)

    schemes = {}
    with open(args.color_schemes) as f:
        for counter, line in enumerate(f, start=1):
            schemes[counter] = line.strip().split("\t")

    with open(args.output, 'w') as f:
        for trait_name, trait_array in assignment.items():
            color_array = schemes[len(trait_array)]
            extra_trait_values = forced_colors.get(trait_name, [])
            extra_color_values = forced_colors.get(trait_name, {}).values()

            zipped = zip(chain(trait_array, extra_trait_values),
                         chain(color_array, extra_color_values))
            for trait_value, color in zipped:
                print("\t".join((trait_name, trait_value, color)), file=f)
            print(file=f)
