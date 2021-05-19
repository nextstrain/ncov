#!/usr/bin/env python3
import argparse
from augur.io import open_file
from augur.utils import read_strains
from collections import defaultdict
import numpy as np
import pandas as pd
import sys
import time


def jaccard_distance(set_a, set_b, missing_from_a):
    """Calculate the Jaccard similarity between two sets of categorical values.
    """
    return 1.0 - (len(set_a & set_b) / (len(set_a | set_b) + missing_from_a))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, help="GISAID metadata with 'aa_substitutions' field for all records")
    parser.add_argument("--focal", required=True, help="strains in the focal data set to find closest strains from the full set")
    parser.add_argument("--substitutions-field", default="aa_substitutions", help="field with annotations of substitutions to use for distance calculations")
    parser.add_argument("--drop-duplicates", action="store_true", help="drop duplicate combinations of substitutions, limiting distances to the first strain with a given set of substitutions in the metadata.")
    parser.add_argument("--output", required=True, help="strain names of closest strains to focal set")

    args = parser.parse_args()

    # Load names of focal strains.
    focal_strains = sorted(read_strains(args.focal))

    # Load full metadata.
    try:
        print("Reading metadata")
        metadata = pd.read_csv(
            args.metadata,
            sep="\t",
            dtype={
                "strain": "string",
                "name": "string",
            },
            usecols=["strain", args.substitutions_field]
        )
    except ValueError as error:
        if args.substitutions_field in str(error):
            error_message = f"ERROR: The substitutions field '{args.substitutions_field}' is not a column in '{args.metadata}'."
        else:
            error_message = f"ERROR: {error}"

        print(error_message, file=sys.stderr)
        sys.exit(1)

    # Drop duplicate mutations.
    if args.drop_duplicates:
        metadata = metadata.drop_duplicates(args.substitutions_field)

    all_strains = metadata["strain"].values
    focal_strains = [strain for strain in focal_strains if strain in all_strains]
    other_strains = [strain for strain in all_strains if strain not in focal_strains]

    print("Parsing and enumerating mutations per strain")
    mutations_by_strain = {}
    index_by_mutation = {}
    current_index = 0
    for strain, value in zip(all_strains, metadata[args.substitutions_field].values):
        strain_mutations = set(str(value)[1:-1].split(","))
        mutations_by_strain[strain] = strain_mutations

        if strain in focal_strains:
            for mutation in strain_mutations:
                if mutation not in index_by_mutation:
                    index_by_mutation[mutation] = current_index
                    current_index += 1

    del metadata

    total_mutations = len(index_by_mutation)
    print(f"Found {total_mutations} distinct mutations in focal strains")

    print("Mapping enumerated mutations to strains")
    mutation_indices_by_strain = defaultdict(set)
    missing_mutations_by_strain = defaultdict(int)
    for strain in all_strains:
        for mutation in mutations_by_strain[strain]:
            if mutation in index_by_mutation:
                mutation_indices_by_strain[strain].add(index_by_mutation[mutation])
            else:
                missing_mutations_by_strain[strain] += 1

    del mutations_by_strain
    del index_by_mutation

    print(f"Calculating distances for {len(other_strains)} non-focal strains")
    with open_file(args.output, "w") as oh:
        oh.write("strain\tclosest strain\tdistance\n")

        then = time.time()
        counter = 0
        for strain in other_strains:
            if counter > 0 and counter % 10000 == 0:
                now = time.time()
                print(f"Processed {counter} records in {now - then:.0f} seconds")
                then = now

            counter += 1
            min_distance = np.inf
            min_strain = None

            for focal_strain in focal_strains:
                distance = jaccard_distance(
                    mutation_indices_by_strain[strain],
                    mutation_indices_by_strain[focal_strain],
                    missing_mutations_by_strain[strain]
                )

                if distance < min_distance:
                    min_distance = distance
                    min_strain = focal_strain

                if min_distance == 0:
                    break

            oh.write("{strain}\t{closest_strain}\t{distance:.4f}\n".format(
                strain=strain,
                closest_strain=min_strain,
                distance=min_distance
            ))
