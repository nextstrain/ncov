import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="input ordering file")
    parser.add_argument('--color-schemes', type=str, required=True, help="input color schemes file")
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
            zipped = list(zip(trait_array, color_array))
            for trait_value, color in zipped:
                f.write(trait_name + "\t" + trait_value + "\t" + color + "\n")
            f.write("\n")
