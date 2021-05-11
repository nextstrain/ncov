import argparse
from augur.utils import read_metadata
import pandas as pd
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="metadata to be sanitized")
    parser.add_argument("--strip-prefixes", nargs="+", help="prefixes to strip from strain names in the metadata")
    parser.add_argument("--output", required=True, help="sanitized metadata")

    args = parser.parse_args()

    metadata, columns = read_metadata(args.metadata)
    metadata = pd.DataFrame.from_dict(metadata, orient="index")

    if args.strip_prefixes:
        prefixes = "|".join(args.strip_prefixes)
        pattern = f"^({prefixes})"

        metadata["strain"] = metadata["strain"].apply(
            lambda strain: re.sub(
                pattern,
                "",
                strain
            )
        )

    metadata.to_csv(
        args.output,
        sep="\t",
        index=False
    )
