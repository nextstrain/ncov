import argparse
from augur.utils import read_metadata
from Bio import SeqIO
import csv

EMPTY = ''

# This script was written in preparation for a future augur where commands
# may take multiple metadata files, thus making this script unnecessary!
#
# Merging logic:
# - Order of supplied TSVs matters
# - All columns are included (i.e. union of all columns present)
# - The last non-empty value read (from different TSVs) is used. I.e. values are overwritten.
# - Missing data is represented by an empty string

def parse_args():
    parser = argparse.ArgumentParser(
        description="Custom script to combine metadata files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', required=True, nargs='+', metavar="TSV", help="Metadata files")
    parser.add_argument('--origins', required=True, nargs='+', metavar="STR", help="Names of origins (metadata columns will be created from these)")
    parser.add_argument('--sequences', required=True, nargs='+', metavar="FASTA", help="Sequences files per-origin")
    parser.add_argument('--output', required=True, metavar="TSV", help="output (merged) metadata")
    args = parser.parse_args()
    return args

def get_sequences_per_origin(origins, fastas):
    assert(len(origins)==len(fastas))
    contents = {}
    for idx in range(0, len(origins)):
        names = set()
        with open(fastas[idx]) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                names.add(record.id)
        contents[origins[idx]] = names
    return contents

if __name__ == '__main__':
    args = parse_args()

    metadata = []
    for fname in args.metadata:
        data, columns = read_metadata(fname)
        metadata.append({'name': fname, 'data': data, 'columns': columns})

    # PARSE FASTA FILES (FOR EACH ORIGIN)
    fasta_contents = get_sequences_per_origin(args.origins, args.sequences)

    # SUMMARISE INPUT METADATA
    print(f"Parsed {len(metadata)} metadata TSVs")
    for m in metadata:
        print(f"\t{m['name']}: {len(m['data'].keys())} strains x {len(m['columns'])} columns")
    print(f"Parsed {len(fasta_contents.keys())} FASTAs")
    for origin, contents in fasta_contents.items():
        print(f"\torigin {origin}: {len(contents)} sequences")

    # BUILD UP COLUMN NAMES TO PRESERVE ORDER
    combined_columns = []
    for m in metadata:
        combined_columns.extend([c for c in m['columns'] if c not in combined_columns])
    combined_columns.extend(list(fasta_contents.keys()))

    # ADD IN VALUES ONE BY ONE, OVERWRITING AS NECESSARY
    combined_data = metadata[0]['data']
    for strain in combined_data:
        for column in combined_columns:
            if column not in combined_data[strain]:
                combined_data[strain][column] = EMPTY

    for m in metadata[1:]:
        for strain, row in m['data'].items():
            if strain not in combined_data:
                combined_data[strain] = {c:EMPTY for c in combined_columns}
            for column in combined_columns:
                if column in row:
                    existing_value = combined_data[strain][column]
                    new_value = row[column]
                    # overwrite _ANY_ existing value if the overwriting value is non empty (and different)!
                    if new_value != EMPTY and new_value != existing_value:
                        if existing_value != EMPTY:
                            print(f"[{strain}::{column}] Overwriting {combined_data[strain][column]} with {new_value}")
                        combined_data[strain][column] = new_value

    ## For each sequence (in metadata) add in boolean value indicating if it's present in each sequence (fasta) file
    orphan_metadata_count = 0
    for strain, data in combined_data.items():
        for origin, names in fasta_contents.items():
            data[origin] = "yes" if strain in names else "no"
        if not any(data[origin] for origin in fasta_contents.keys()):
            orphan_metadata_count+=1        
    

    print(f"Combined metadata: {len(combined_data.keys())} strains x {len(combined_columns)} columns")
    if orphan_metadata_count:
        print(f"\tWARNING: {orphan_metadata_count} samples were missing from supplied sequence files")

    with open(args.output, 'w') as fh:
        tsv_writer = csv.writer(fh, delimiter='\t')
        tsv_writer.writerow(combined_columns)
        for row in combined_data.values():
            tsv_writer.writerow([row[column] for column in combined_columns])