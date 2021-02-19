import argparse
from augur.utils import read_metadata
from Bio import SeqIO
import csv
import sys

EMPTY = ''

# This script was written in preparation for a future augur where commands
# may take multiple metadata files, thus making this script unnecessary!
#
# Merging logic:
# - Order of supplied TSVs matters
# - All columns are included (i.e. union of all columns present)
# - The last non-empty value read (from different TSVs) is used. I.e. values are overwritten.
# - Missing data is represented by an empty string
#
# We use one-hot encoding to specify which origin(s) a piece of metadata came from

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Custom script to combine metadata files from different origins.
        In the case where metadata files specify different values, the latter provided file will take priority.
        Columns will be added for each origin with values "yes" or "no" to identify the input source (origin) of each sample.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', required=True, nargs='+', metavar="TSV", help="Metadata files")
    parser.add_argument('--origins', required=True, nargs='+', metavar="STR", help="Names of origins (order should match provided metadata)")
    parser.add_argument('--output', required=True, metavar="TSV", help="Output (merged) metadata")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    try:
        assert(len(args.metadata)==len(args.origins))
        assert(len(args.origins)>1)
    except AssertionError:
        print("Error. Please check your inputs - there must be the same number of metadata files as origins provided, and there must be more than one of each!")
        sys.exit(2)

    # READ IN METADATA FILES
    metadata = []
    for (origin, fname) in zip(args.origins, args.metadata):
        data, columns = read_metadata(fname)
        metadata.append({'origin': origin, "fname": fname, 'data': data, 'columns': columns, 'strains': {s for s in data.keys()}})

    # SUMMARISE INPUT METADATA
    print(f"Parsed {len(metadata)} metadata TSVs")
    for m in metadata:
        print(f"\t{m['origin']} ({m['fname']}): {len(m['data'].keys())} strains x {len(m['columns'])} columns")

    # BUILD UP COLUMN NAMES FROM MULTIPLE INPUTS TO PRESERVE ORDER
    combined_columns = []
    for m in metadata:
        combined_columns.extend([c for c in m['columns'] if c not in combined_columns])
    combined_columns.extend(list(args.origins))

    # ADD IN VALUES ONE BY ONE, OVERWRITING AS NECESSARY
    combined_data = metadata[0]['data']
    for strain in combined_data:
        for column in combined_columns:
            if column not in combined_data[strain]:
                combined_data[strain][column] = EMPTY    
    
    for idx in range(1, len(metadata)):
        for strain, row in metadata[idx]['data'].items():
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

    # one-hot encoding for origin
    # note that we use "yes" / "no" here as Booleans are problematic for `augur filter`
    for metadata_entry in metadata:
        origin = metadata_entry['origin']
        for strain in combined_data:
            combined_data[strain][origin] = "yes" if strain in metadata_entry['strains'] else "no"

    print(f"Combined metadata: {len(combined_data.keys())} strains x {len(combined_columns)} columns")

    with open(args.output, 'w') as fh:
        tsv_writer = csv.writer(fh, delimiter='\t')
        tsv_writer.writerow(combined_columns)
        for row in combined_data.values():
            tsv_writer.writerow([row[column] for column in combined_columns])