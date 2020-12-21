"""Split sequences for multiple sequence alignment into smaller chunks to speed up alignment.
"""
import argparse
from augur.align import read_sequences
from Bio import SeqIO
from pathlib import Path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences to partition into smaller chunks")
    parser.add_argument("--sequences-per-group", required=True, type=int, help="number of sequences to include in each group")
    parser.add_argument("--output-dir", required=True, help="directory to write out partitioned sequences")

    args = parser.parse_args()

    # Read sequences with augur to benefit from additional checks for duplicates.
    sequences = SeqIO.parse(args.sequences, "fasta")

    # Create the requested output directory.
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    i = 0
    group_number = 0
    for sequence in sequences:
        if i == args.sequences_per_group:
            print(f"Write out {i} sequences for partition {group_number}")
            out_handle.close()
            i = 0
            group_number += 1

        if i == 0:
            output_path = str(output_dir / Path(f"{group_number}.fasta"))
            out_handle = open(output_path, "w")

        SeqIO.write(sequence, out_handle, "fasta")
        i += 1
    else:
        print(f"Write out {i} sequences for partition {group_number}")
        out_handle.close()
