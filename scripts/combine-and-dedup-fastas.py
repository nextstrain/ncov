import argparse
from augur.io import read_sequences
from Bio import SeqIO
import hashlib
import sys
import textwrap


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Combine and dedup FASTAs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str,  nargs="+", metavar="FASTA", required=True, help="input FASTAs")
    parser.add_argument('--warn-about-duplicates', action="store_true", help="warn the user about duplicate strains instead of exiting with an error. The output will include the first occurrence of a duplicate sequence.")
    parser.add_argument('--output', type=str, metavar="FASTA", required=True, help="output FASTA")
    args = parser.parse_args()

    sequence_hash_by_name = {}
    duplicate_strains = set()

    counter = 0
    with open(args.output, "w") as output_handle:
        # Stream sequences from all input files into a single output file,
        # skipping duplicate records (same strain and sequence) and noting
        # mismatched sequences for the same strain name.
        for record in read_sequences(*args.input):
            counter += 1
            if counter % 10000 == 0:
                print(f"Processed {counter} records")

            # Hash each sequence and check whether another sequence with the
            # same name already exists and if the hash is different.
            sequence_hash = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
            if record.name in sequence_hash_by_name:
                # If the hashes differ (multiple entries with the same
                # strain name but different sequences), we keep the first
                # sequence and add the strain to a list of duplicates to
                # report at the end.
                if sequence_hash_by_name.get(record.name) != sequence_hash:
                    duplicate_strains.add(record.name)

                # If the current strain has been seen before, don't write
                # out its sequence again.
                continue

            sequence_hash_by_name[record.name] = sequence_hash
            SeqIO.write(record, output_handle, 'fasta')

    if len(duplicate_strains) > 0:
        error_mode = "ERROR"
        exit_code = 1

        if args.warn_about_duplicates:
            error_mode = "WARNING"
            exit_code = 0

        print(
            f"{error_mode}: Detected the following duplicate input strains with different sequences:",
            file=sys.stderr
        )
        for strain in duplicate_strains:
            print(textwrap.indent(strain, "    "), file=sys.stderr)

        if not args.warn_about_duplicates:
            print(
                "Use the `--warn-about-duplicates` flag, if you prefer to accept",
                "the first occurrence of duplicate sequences based on the order",
                "of the given input sequences.",
                file=sys.stderr
            )

        sys.exit(exit_code)
