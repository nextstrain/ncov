"""Split sequences for multiple sequence alignment into smaller chunks to speed up alignment.
"""
import argparse
from augur.align import read_sequences
from Bio import SeqIO
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, nargs="+", help="FASTA file of sequences to partition into smaller chunks")
    parser.add_argument("--sequences-per-group", required=True, type=int, help="number of sequences to include in each group")
    parser.add_argument("--output-dir", required=True, help="directory to write out partitioned sequences")

    args = parser.parse_args()

    # Read sequences with augur to benefit from additional checks for duplicates.
    sequences = list(read_sequences(*args.sequences).values())

    # Create the requested output directory.
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Determine partition indices.
    indices = list(range(0, len(sequences), args.sequences_per_group))

    # Add a final index to represent the end of the last partition.
    if indices[-1] != len(sequences):
        indices.append(len(sequences))

    # Partition sequences into groups of no more than the requested number.
    for i in range(len(indices) - 1):
        # Save partitioned sequences to a new FASTA file named after the partition number.
        print("Write out %i sequences for partition %i" % (indices[i + 1] - indices[i], i))
        output_path = output_dir / Path("%i.fasta" % i)
        SeqIO.write(sequences[indices[i]:indices[i + 1]], str(output_path), 'fasta')
