from Bio import SeqIO
import argparse
from augur.align import read_sequences

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Combine and dedup FASTAs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str,  nargs="+", metavar="FASTA", required=True, help="input FASTAs")
    parser.add_argument('--output', type=str, metavar="FASTA", required=True, help="output FASTA")
    args = parser.parse_args()

    # Read sequences with augur to benefit from additional checks for duplicates.
    sequences = read_sequences(*args.input)

    # Convert dictionary of sequences by id to a list, for compatibility with
    # augur versions <9.0.0.
    if isinstance(sequences, dict):
        sequences = list(sequences.values())

    SeqIO.write(sequences, args.output, 'fasta')
