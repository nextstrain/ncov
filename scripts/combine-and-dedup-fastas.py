import argparse
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
    parser.add_argument('--output', type=str, metavar="FASTA", required=True, help="output FASTA")
    args = parser.parse_args()

    sequence_hash_by_name = {}
    duplicate_strains = set()

    counter = 0
    with open(args.output, "w") as output_handle:
        for filename in args.input:
            for record in SeqIO.parse(filename, 'fasta'):
                counter += 1
                if counter % 10000 == 0:
                    print(f"Processed {counter} records")

                sequence_hash = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
                if record.name in sequence_hash_by_name and sequence_hash_by_name.get(record.name) != sequence_hash:
                    duplicate_strains.add(record.name)
                    continue

                sequence_hash_by_name[record.name] = sequence_hash
                SeqIO.write(record, output_handle, 'fasta')

    if len(duplicate_strains) > 0:
        print(
            "WARNING: Detected the following duplicate input strains with different sequences:",
            file=sys.stderr
        )
        for strain in duplicate_strains:
            print(textwrap.indent(strain, "    "), file=sys.stderr)
