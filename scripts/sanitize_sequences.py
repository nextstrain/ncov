import argparse
from augur.io import open_file, read_sequences, write_sequences
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sequences", nargs="+", required=True, help="sequences to be sanitized")
    parser.add_argument("--strip-prefixes", nargs="+", help="prefixes to strip from strain names in the sequences")
    parser.add_argument("--output", required=True, help="sanitized sequences")

    args = parser.parse_args()

    if args.strip_prefixes:
        prefixes = "|".join(args.strip_prefixes)
        pattern = f"^({prefixes})"
    else:
        pattern = ""

    with open_file(args.output, "w") as output_handle:
        for sequence in read_sequences(*args.sequences):
            sequence.id = re.sub(pattern, "", sequence.id)
            write_sequences(sequence, output_handle)
