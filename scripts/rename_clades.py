import argparse
import yaml

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rename clades in clades.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-clade-files', type=str, nargs='+', required=True, help="input clade files")
    parser.add_argument('--name-mapping', type=str, required=False, help="YAML mapping between Nextstrain clades and display names")
    parser.add_argument('--output-clades', type=str, required=True, help="renamed clade file")
    args = parser.parse_args()

    # read name mapping from input yaml file
    if args.name_mapping:
        with open(args.name_mapping) as fh:
            name_mapping = yaml.load(fh, Loader=yaml.FullLoader)
    else:
        name_mapping = {}


    # write output into one consolidated file
    out_clades = open(args.output_clades, "w")

    # loop over input file and replace clade names were appropriate line by line
    for fname in args.input_clade_files:
        with open(fname) as fh:
            for line in fh:
                fields = line.strip('\n').split('\t')
                if len(fields) < 3:
                    continue
                fields[0] = name_mapping.get(fields[0], fields[0])
                # if clade definition is based on other clade, replace name
                if fields[1]=='clade':
                    fields[2] = name_mapping.get(fields[2], fields[2])
                out_clades.write('\t'.join(fields)+'\n')

    out_clades.close()
