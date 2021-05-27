import argparse
from pathlib import Path
import pandas as pd
import re
import sys

from utils import extract_tar_file_contents

# Define all possible geographic scales we could expect in the GISAID location
# field.
LOCATION_FIELDS = (
    "region",
    "country",
    "division",
    "location",
)

# Define possible strain name fields.
STRAIN_FIELDS = (
    "strain",
    "name",
    "Virus name",
)

ACCESSION_FIELDS = (
    "gisaid_epi_isl",
    "genbank_accession",
)


def parse_location_string(location_string, location_fields):
    """Parse location string from GISAID into the given separate geographic scales
    and return a dictionary of parse values by scale.

    Parameters
    ----------
    location_string : str
    location_fields : list

    Returns
    -------
    dict :
        dictionary of geographic fields parsed from the given string

    >>> location_fields = ["region", "country", "division", "location"]
    >>> parse_location_string("Asia / Japan", location_fields)
    {'region': 'Asia', 'country': 'Japan', 'division': '?', 'location': '?'}

    >>> parse_location_string("Europe / Iceland / Reykjavik", location_fields)
    {'region': 'Europe', 'country': 'Iceland', 'division': 'Reykjavik', 'location': '?'}

    >>> parse_location_string("North America / USA / Washington / King County", location_fields)
    {'region': 'North America', 'country': 'USA', 'division': 'Washington', 'location': 'King County'}

    Additional location entries beyond what has been specified should be stripped from output.

    >>> parse_location_string("North America / USA / Washington / King County / Extra field", location_fields)
    {'region': 'North America', 'country': 'USA', 'division': 'Washington', 'location': 'King County'}

    Trailing location delimiters should be stripped from the output.

    >>> parse_location_string("North America / USA / Washington / King County / ", location_fields)
    {'region': 'North America', 'country': 'USA', 'division': 'Washington', 'location': 'King County'}

    Handle inconsistently delimited strings.

    >>> parse_location_string("North America/USA/New York/New York", location_fields)
    {'region': 'North America', 'country': 'USA', 'division': 'New York', 'location': 'New York'}
    >>> parse_location_string("Europe/ Lithuania", location_fields)
    {'region': 'Europe', 'country': 'Lithuania', 'division': '?', 'location': '?'}

    """
    # Try to extract values for specific geographic scales.
    values = re.split(r"[ ]*/[ ]*", location_string)

    # Create a default mapping of location fields to missing values and update
    # these from the values in the location string.
    locations = {field: "?" for field in location_fields}
    locations.update(dict(zip(location_fields, values)))

    return locations


def resolve_duplicates(metadata, strain_field, error_on_duplicates=False):
    """Resolve duplicate records for a given strain field and return a deduplicated
    data frame. This approach chooses the record with the most recent database
    accession, if accession fields exist, or the first record for a given strain
    name. Optionally, raises an error when duplicates are detected, reporting
    the list of those duplicate records.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Metadata table that may or may not have duplicate records by the given strain field.

    strain_field : string
        Name of the metadata field corresponding to the strain name and which is used to identify duplicates.

    error_on_duplicates : boolean
        Whether to throw an error on detection of duplicates or not.

    Returns
    -------
    pandas.DataFrame :
        Metadata with no duplicate records by the given strain field.

    Raises
    ------
    ValueError
        If duplicates are detected and `error_on_duplicates` is `True`.
    """
    # Create a list of duplicate records by strain name.
    duplicates = metadata.loc[
        metadata.duplicated(strain_field),
        strain_field
    ].values

    if len(duplicates) == 0:
        # No duplicates, so return the original metadata.
        return metadata

    if error_on_duplicates:
        # Duplicates and error requested on duplicates, so throw an exception
        # with the list of duplicate strains.
        raise ValueError(", ".join(duplicates))

    # Try to resolve the duplicates by preferring records with the most recent
    # database accession. First, check for standard accession fields.
    accession_fields = [
        field
        for field in ACCESSION_FIELDS
        if field in metadata.columns
    ]

    # If any of these fields exists, sort by strain name and accessions in
    # ascending order and take the last record (the most recent accession for a
    # given strain). Otherwise, sort and group by strain name and take the last
    # record from each group. It is possible for metadata to contain fields for
    # multiple database accessions and for these fields to be incomplete for
    # some databases (for example, GISAID accession is `?` but GenBank accession
    # is not). By sorting across all fields, we use information from all
    # available accession fields. If all fields contain missing values (e.g.,
    # "?"), we end up returning the last record for a given strain as a
    # reasonable default.
    sort_fields = [strain_field]
    if len(accession_fields) > 0:
        sort_fields.extend(accession_fields)

    # Return the last record from each group after sorting by strain and
    # available accessions.
    return metadata.sort_values(sort_fields).drop_duplicates(strain_field, "last")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Sanitize metadata from different sources, applying operations in the same order they appear in the full help (-h).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="metadata to be sanitized")
    parser.add_argument("--parse-location-field", help="split the given GISAID location field on '/' and create new columns for region, country, etc. based on available data. Replaces missing geographic data with '?' values.")
    parser.add_argument("--rename-fields", nargs="+", help="rename specific fields from the string on the left of the equal sign to the string on the right (e.g., 'Virus name=strain')")
    parser.add_argument("--strip-prefixes", nargs="+", help="prefixes to strip from strain names in the metadata")
    parser.add_argument("--error-on-duplicate-strains", action="store_true", help="exit with an error if any duplicate strains are found. By default, duplicates are resolved by preferring most recent accession id or the last record.")
    parser.add_argument("--output", required=True, help="sanitized metadata")

    args = parser.parse_args()

    # If the input is a tarball, try to find a metadata file inside the archive.
    metadata_file = args.metadata
    tar_handle = None
    if ".tar" in Path(args.metadata).suffixes:
        try:
            metadata_file, tar_handle = extract_tar_file_contents(
                args.metadata,
                "metadata"
            )
        except FileNotFoundError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            sys.exit(1)

    # Read metadata with pandas because Augur's read_metadata utility does not
    # support metadata without a "strain" or "name" field.
    metadata = pd.read_csv(
        metadata_file,
        sep=None,
        engine="python",
        skipinitialspace=True,
        dtype={
            "strain": "string",
            "name": "string",
        }
    ).fillna("")

    # Close tarball after reading metadata if it is open still.
    if tar_handle is not None:
        tar_handle.close()

    if args.parse_location_field and args.parse_location_field in metadata.columns:
        # Parse GISAID location field into separate fields for geographic
        # scales. Replace missing field values with "?".
        locations = pd.DataFrame(
            (
                parse_location_string(location, LOCATION_FIELDS)
                for location in metadata[args.parse_location_field].values
            )
        )

        # Combine new location columns with original metadata and drop the
        # original location column.
        metadata = pd.concat(
            [
                metadata,
                locations
            ],
            axis=1
        ).drop(columns=[args.parse_location_field])

    new_column_names = {}
    if args.rename_fields:
        # Rename specific columns using rules like "Virus name=strain".
        for rule in args.rename_fields:
            if "=" in rule:
                old_column, new_column = rule.split("=")
                new_column_names[old_column] = new_column
            else:
                print(
                    f"WARNING: missing mapping of old to new column in form of 'Virus name=strain' for rule: '{rule}'.",
                    file=sys.stderr
                )

    # Rename columns as needed.
    if len(new_column_names) > 0:
        metadata = metadata.rename(columns=new_column_names)

    # Determine field for strain name.
    strain_field = None
    for field in STRAIN_FIELDS:
        if field in metadata.columns:
            strain_field = field
            break

    if strain_field is None:
        print(
            f"ERROR: None of the available columns match possible strain name fields ({', '.join(STRAIN_FIELDS)}).",
            f"Available columns are: {metadata.columns.values}",
            file=sys.stderr
        )
        sys.exit(1)

    if args.strip_prefixes:
        prefixes = "|".join(args.strip_prefixes)
        pattern = f"^({prefixes})"

        metadata[strain_field] = metadata[strain_field].apply(
            lambda strain: re.sub(
                pattern,
                "",
                strain
            )
        )

    # Remove whitespaces from strain names since they are not allowed in FASTA
    # record names.
    metadata[strain_field] = metadata[strain_field].str.replace(" ", "")

    # Check for duplicates and try to resolve these by default.
    try:
        metadata = resolve_duplicates(
            metadata,
            strain_field,
            error_on_duplicates=args.error_on_duplicate_strains
        )
    except ValueError as e:
        print(f"ERROR: The following strains have duplicate records: {e}", file=sys.stderr)
        sys.exit(1)

    metadata.to_csv(
        args.output,
        sep="\t",
        index=False
    )
