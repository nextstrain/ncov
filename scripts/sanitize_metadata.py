import argparse
from augur.io import open_file, read_metadata
import csv
import os
from pathlib import Path
import pandas as pd
import re
import shutil
import sys
from tempfile import NamedTemporaryFile

from utils import extract_tar_file_contents

# Define all possible geographic scales we could expect in the GISAID location
# field.
LOCATION_FIELDS = (
    "region",
    "country",
    "division",
    "location",
)


class MissingColumnException(Exception):
    """An exception caused by a missing column that was expected in the metadata.

    """
    pass


class DuplicateException(Exception):
    """An exception caused by the presence of any duplicate metadata records by
    strain name.

    """
    pass


def parse_new_column_names(renaming_rules):
    """Parse the mapping of current to new column names from the given list of renaming rules.

    Parameters
    ----------
    renaming_rules : list[str]
        A list of strings mapping an old column name to a new one delimited by an equal symbol (e.g., "old_column=new_column").

    Returns
    -------
    dict :
        A mapping of new column names for each old column name.

    >>> parse_new_column_names(["old=new", "new=old"])
    {'old': 'new', 'new': 'old'}
    >>> parse_new_column_names(["old->new"])
    {}

    """
    new_column_names = {}
    for rule in renaming_rules:
        if "=" in rule:
            old_column, new_column = rule.split("=")
            new_column_names[old_column] = new_column
        else:
            print(
                f"WARNING: missing mapping of old to new column in form of 'Virus name=strain' for rule: '{rule}'.",
                file=sys.stderr
            )

    return new_column_names


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


def strip_prefixes(strain_name, prefixes):
    """Strip the given prefixes from the given strain name.

    Parameters
    ----------
    strain_name : str
        Name of a strain to be sanitized
    prefixes : list[str]
        A list of prefixes to be stripped from the strain name.

    Returns
    -------
    str :
        Strain name without any of the given prefixes.


    >>> strip_prefixes("hCoV-19/RandomStrain/1/2020", ["hCoV-19/", "SARS-CoV-2/"])
    'RandomStrain/1/2020'
    >>> strip_prefixes("SARS-CoV-2/RandomStrain/2/2020", ["hCoV-19/", "SARS-CoV-2/"])
    'RandomStrain/2/2020'
    >>> strip_prefixes("hCoV-19/RandomStrain/1/2020", ["SARS-CoV-2/"])
    'hCoV-19/RandomStrain/1/2020'

    """
    joined_prefixes = "|".join(prefixes)
    pattern = f"^({joined_prefixes})"
    return re.sub(pattern, "", strain_name)


def get_database_ids_by_strain(metadata_file, metadata_id_columns, database_id_columns, metadata_chunk_size, error_on_duplicates=False):
    """Get a mapping of all database ids for each strain name.

    Parameters
    ----------
    metadata_file : str or Path-like or file object
        Path or file object for a metadata file to process.
    metadata_id_columns : list[str]
        A list of potential id columns for strain names in the metadata.
    database_id_columns : list[str]
        A list of potential database id columns whose values can be used to deduplicate records with the same strain name.
    metadata_chunk_size : int
        Number of records to read into memory at once from the metadata.
    error_on_duplicates : bool
        Throw an error when duplicate records are detected.

    Returns
    -------
    str or Path-like or file object or None :
        Path or file object containing the mapping of database ids for each
        strain name (one row per combination). Returns None, if no valid
        database ids were found and no duplicates exist.

    Raises
    ------
    DuplicateException :
        When duplicates are detected and the caller has requested an error on duplicates.
    MissingColumnException :
        When none of the requested metadata id columns exist.

    """
    try:
        metadata_reader = read_metadata(
            metadata_file,
            id_columns=metadata_id_columns,
            chunk_size=metadata_chunk_size,
        )
    except Exception as error:
        # Augur's `read_metadata` function can throw a generic Exception when
        # the input is missing id columns. This exception is not easily
        # distinguished from any other error, so we check the contents of the
        # error message and raise a more specific error for better handling of
        # unexpected errors.
        if "None of the possible id columns" in str(error):
            raise MissingColumnException(str(error)) from error
        else:
            raise

    # Track strains we have observed, so we can alert the caller to duplicate
    # strains when an error on duplicates has been requested.
    observed_strains = set()
    duplicate_strains = set()

    with NamedTemporaryFile(delete=False, mode="wt", encoding="utf-8", newline="") as mapping_file:
        mapping_path = mapping_file.name
        header = True

        for metadata in metadata_reader:
            # Check for database id columns.
            valid_database_id_columns = metadata.columns.intersection(
                database_id_columns
            )
            if mapping_path and len(valid_database_id_columns) == 0:
                # Do not write out mapping of ids. Default to error on
                # duplicates, since we have no way to resolve duplicates
                # automatically.
                mapping_path = None
                error_on_duplicates = True
                print(
                    "WARNING: Skipping deduplication of metadata records.",
                    f"None of the possible database id columns ({database_id_columns}) were found in the metadata's columns {tuple([metadata.index.name] + metadata.columns.values.tolist())}",
                    file=sys.stderr
                )

            # Track duplicates in memory, as needed.
            if error_on_duplicates:
                for strain in metadata.index.values:
                    if strain in observed_strains:
                        duplicate_strains.add(strain)
                    else:
                        observed_strains.add(strain)

            if mapping_path:
                # Write mapping of database and strain ids to disk.
                metadata.loc[:, valid_database_id_columns].to_csv(
                    mapping_file,
                    sep="\t",
                    header=header,
                    index=True,
                )
                header = False

    # Clean up temporary file.
    if mapping_path is None:
        os.unlink(mapping_file.name)

    if error_on_duplicates and len(duplicate_strains) > 0:
        duplicates_file = metadata_file + ".duplicates.txt"
        with open(duplicates_file, "w") as oh:
            for strain in duplicate_strains:
                oh.write(f"{strain}\n")

        raise DuplicateException(f"{len(duplicate_strains)} strains have duplicate records. See '{duplicates_file}' for more details.")

    return mapping_path


def filter_duplicates(metadata, database_ids_by_strain):
    """Filter duplicate records by the strain name in the given data frame index
    using the given file containing a mapping of strain names to database ids.

    Database ids allow us to identify duplicate records that need to be
    excluded. We prefer the latest record for a given strain name across all
    possible database ids and filter out all other records for that same strain
    name.

    Parameters
    ----------
    metadata : pandas.DataFrame
        A data frame indexed by strain name.
    database_ids_by_strain : str or Path-like or file object
        Path or file object containing the mapping of database ids for each strain name (one row per combination).

    Returns
    -------
    pandas.DataFrame :
        A filtered data frame with no duplicate records.

    """
    # Get strain names for the given metadata.
    strain_ids = set(metadata.index.values)

    # Get the mappings of database ids to strain names for the current strains.
    with open(database_ids_by_strain, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        # The mapping file stores the strain name in the first column. All other
        # fields are database ids.
        strain_field = reader.fieldnames[0]
        database_id_columns = reader.fieldnames[1:]

        # Keep only records matching the current strain ids.
        mappings = pd.DataFrame([
            row
            for row in reader
            if row[strain_field] in strain_ids
        ])

    # Check for duplicate strains in the given metadata or strains that do not
    # have any mappings. If there are none, return the metadata as it is. If
    # duplicates or strains without mappings exist, filter them out.
    if any(mappings.duplicated(strain_field)) or len(strain_ids) != mappings.shape[0]:
        # Create a list of database ids of records to keep. To this end, we sort by
        # database ids in descending order such that the latest record appears
        # first, then we take the first record for each strain name.
        records_to_keep = mappings.sort_values(
            database_id_columns,
            ascending=False
        ).groupby(strain_field).first()

        # Select metadata corresponding to database ids to keep. Database ids
        # may not be unique for different strains (e.g., "?"), so we need to
        # merge on strain name and database ids. Additionally, the same strain
        # may appear multiple times in the metadata with the same id. These
        # accidental duplicates will also produce a merge with records to keep
        # that is not a one-to-one merge. To handle this case, we need to drop
        # any remaining duplicate records by strain name. The order that we
        # resolve these duplicates does not matter, since the fields we would
        # use to resolve these contain identical values.
        merge_columns = sorted(set([strain_field]) | set(database_id_columns))
        metadata = metadata.reset_index().merge(
            records_to_keep,
            on=merge_columns,
        ).drop_duplicates(subset=strain_field).set_index(strain_field)

    # Track strains that we've processed and drop these from the mappings file.
    # In this way, we can track strains that have been processed across multiple
    # chunks of metadata and avoid emiting duplicates that appear in different
    # chunks.
    with open(database_ids_by_strain, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        with NamedTemporaryFile(delete=False, mode="wt", encoding="utf-8", newline="") as new_mapping_file:
            new_mapping_path = new_mapping_file.name
            writer = csv.DictWriter(
                new_mapping_file,
                fieldnames=reader.fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in reader:
                if row[strain_field] not in strain_ids:
                    writer.writerow(row)

    # After writing out the new mapping of ids without strains we just
    # processed, copy the new mapping over the original file and delete the
    # temporary new mapping file.
    shutil.copyfile(
        new_mapping_path,
        database_ids_by_strain,
    )
    os.unlink(new_mapping_path)

    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage="Sanitize metadata from different sources, applying operations (deduplicate, parse location field, strip prefixes, and rename fields) in the same order they appear in the full help (-h).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="metadata to be sanitized")
    parser.add_argument("--metadata-id-columns", default=["strain", "name", "Virus name"], nargs="+", help="names of valid metadata columns containing identifier information like 'strain' or 'name'")
    parser.add_argument("--database-id-columns", default=["Accession ID", "gisaid_epi_isl", "genbank_accession"], nargs="+", help="names of metadata columns that store external database ids for each record (e.g., GISAID, GenBank, etc.) that can be used to deduplicate metadata records with the same strain names.")
    parser.add_argument("--metadata-chunk-size", type=int, default=100000, help="maximum number of metadata records to read into memory at a time. Increasing this number can speed up filtering at the cost of more memory used.")
    parser.add_argument("--error-on-duplicate-strains", action="store_true", help="exit with an error if any duplicate strains are found. By default, duplicates are resolved by preferring most recent accession id or the last record.")
    parser.add_argument("--parse-location-field", help="split the given GISAID location field on '/' and create new columns for region, country, etc. based on available data. Replaces missing geographic data with '?' values.")
    parser.add_argument("--strip-prefixes", nargs="+", help="prefixes to strip from strain names in the metadata")
    parser.add_argument("--rename-fields", nargs="+", help="rename specific fields from the string on the left of the equal sign to the string on the right (e.g., 'Virus name=strain')")
    parser.add_argument("--output", required=True, help="sanitized metadata")

    args = parser.parse_args()

    # Get user-defined metadata id columns to look for.
    metadata_id_columns = args.metadata_id_columns

    # Get user-defined database id columns to use for deduplication.
    database_id_columns = args.database_id_columns

    # If the input is a tarball, try to find a metadata file inside the archive.
    metadata_file = args.metadata
    metadata_is_temporary = False
    if ".tar" in Path(args.metadata).suffixes:
        try:
            temporary_dir, metadata_file = extract_tar_file_contents(
                args.metadata,
                "metadata"
            )
            metadata_is_temporary = True
        except FileNotFoundError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            sys.exit(1)

    # In the first pass through the metadata, map strain names to database ids.
    # We will use this mapping to deduplicate records in the second pass.
    # Additionally, this pass checks for missing id columns and the presence of
    # any duplicate records, in case the user has requested an error on
    # duplicates.
    try:
        database_ids_by_strain = get_database_ids_by_strain(
            metadata_file,
            metadata_id_columns,
            database_id_columns,
            args.metadata_chunk_size,
            args.error_on_duplicate_strains,
        )
    except (DuplicateException, MissingColumnException) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)

    # Parse mapping of old column names to new.
    rename_fields = args.rename_fields if args.rename_fields else []
    new_column_names = parse_new_column_names(rename_fields)

    # In the second pass through the metadata, filter duplicate records,
    # transform records with requested sanitizer steps, and stream the output to
    # disk.
    metadata_reader = read_metadata(
        metadata_file,
        id_columns=metadata_id_columns,
        chunk_size=args.metadata_chunk_size,
    )
    emit_header = True

    with open_file(args.output, "w") as output_file_handle:
        for metadata in metadata_reader:
            if database_ids_by_strain:
                # Filter duplicates.
                metadata = filter_duplicates(
                    metadata,
                    database_ids_by_strain,
                )

            # Reset the data frame index, to make the "strain" column available
            # for transformation.
            strain_field = metadata.index.name
            metadata = metadata.reset_index()

            # Parse GISAID location field into separate fields for geographic
            # scales. Replace missing field values with "?".
            if args.parse_location_field and args.parse_location_field in metadata.columns:
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

            # Strip prefixes from strain names.
            if args.strip_prefixes:
                metadata[strain_field] = metadata[strain_field].apply(
                    lambda strain: strip_prefixes(strain, args.strip_prefixes)
                )

            # Replace whitespaces from strain names with underscores to match GISAID's
            # convention since whitespaces are not allowed in FASTA record names.
            metadata[strain_field] = metadata[strain_field].str.replace(" ", "_")

            # Rename columns as needed, after transforming strain names. This
            # allows us to avoid keeping track of a new strain name field
            # provided by the user.
            if len(new_column_names) > 0:
                metadata = metadata.rename(columns=new_column_names)

            # Write filtered and transformed metadata to the output file.
            metadata.to_csv(
                output_file_handle,
                sep="\t",
                index=False,
                header=emit_header,
            )
            emit_header = False

    if database_ids_by_strain:
        # Delete the database/strain id mapping.
        os.unlink(database_ids_by_strain)

    # Clean up temporary directory and files that came from a tarball.
    if metadata_is_temporary:
        print(f"Cleaning up temporary files in {temporary_dir.name}", file=sys.stderr)
        temporary_dir.cleanup()
