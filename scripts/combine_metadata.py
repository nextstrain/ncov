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
    parser.add_argument('--input', required=True, nargs='+', metavar="TSV", help="Metadata files")
    parser.add_argument('--output', required=True, metavar="TSV", help="output (merged) metadata")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()

    if len(args.input)==1:
        import sys
        from shutil import copyfile
        copyfile(args.input[0], args.output)
        sys.exit(0)

    metadata = []
    for fname in args.input:
        data, columns = read_metadata(fname)
        metadata.append({'name': fname, 'data': data, 'columns': columns})
    
    # SUMMARISE
    print(f"Parsed {len(metadata)} metadata TSVs")
    for m in metadata:
        print(f"\t{m['name']}: {len(m['data'].keys())} strains x {len(m['columns'])} columns")

    # BUILD UP COLUMN NAMES TO PRESERVE ORDER
    combined_columns = []
    for m in metadata:
        combined_columns.extend([c for c in m['columns'] if c not in combined_columns])

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

    print(f"Combined metadata: {len(combined_data.keys())} strains x {len(combined_columns)} columns")


    with open(args.output, 'w') as fh:
        tsv_writer = csv.writer(fh, delimiter='\t')
        tsv_writer.writerow(combined_columns)
        for row in combined_data.values():
            tsv_writer.writerow([row[column] for column in combined_columns])

