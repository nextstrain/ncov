"""Annotate a metadata file with the given sequence index.
"""
import argparse
from augur.io import read_metadata
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="metadata to annotate")
    parser.add_argument("--sequence-index", required=True, help="sequence index from augur index")
    parser.add_argument("--output", required=True, help="metadata annotated with sequence index columns including a 'length' column based on the number of A, C, G, and T bases.")

    args = parser.parse_args()

    metadata = read_metadata(args.metadata)

    index = pd.read_csv(
        args.sequence_index,
        sep="\t",
        dtype={"strain": "string"},
        index_col="strain",
    ).drop(
        columns=["length"],
    )
    index["length"] = index.loc[:, ["A", "C", "G", "T"]].sum(axis=1)
    new_columns = {
        column: f"_{column}"
        for column in index.columns
    }
    index = index.rename(columns=new_columns)

    metadata.join(
        index,
    ).to_csv(
        args.output,
        sep="\t",
    )
