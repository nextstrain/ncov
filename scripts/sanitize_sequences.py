import argparse
from augur.io import open_file, read_sequences, write_sequences
import hashlib
from pathlib import Path
import re
import sys

from utils import stream_tar_file_contents


class DuplicateSequenceError(ValueError):
    pass


def rename_sequences(sequences, pattern):
    """Rename the given sequences' ids by replacing the given patterns with the
    empty string.

    """
    for sequence in sequences:
        # Replace the given patterns in the sequence description with the empty
        # string. For a simple FASTA record with only an identifier in the
        # defline, the description is identical to the `id` and `name` fields.
        # For a complex FASTA record that has spaces in the identifier or other
        # additional information, we need to parse the description to get any
        # trailing components of the strain name.
        sequence.id = re.sub(pattern, "", sequence.description)

        # The name field stores the same information for a simple FASTA input, so we need to override its value, too.
        sequence.name = sequence.id

        # Do not keep additional information that follows the sequence identifier.
        sequence.description = ""

        yield sequence


def drop_duplicate_sequences(sequences, error_on_duplicates=False):
    """Identify and drop duplicate sequences from the given iterator.

    Parameters
    ----------
    sequences : Iterator

    Yields
    ------
    Bio.SeqIO.Seq :
        Unique sequence records

    Raises
    ------
    DuplicateSequenceError :
        If `error_on_duplicates` is True and any duplicate records are found,
        raises an exception with a comma-delimited list of duplicates as the
        message.

    """
    sequence_hash_by_name = {}
    duplicate_strains = set()

    for sequence in sequences:
        # Hash each sequence and check whether another sequence with the same
        # name already exists and if the hash is different.
        sequence_hash = hashlib.sha256(str(sequence.seq).encode("utf-8")).hexdigest()
        if sequence.name in sequence_hash_by_name:
            # If the hashes differ (multiple entries with the same strain name
            # but different sequences), we keep the first sequence and add the
            # strain to a list of duplicates to report at the end.
            if sequence_hash_by_name.get(sequence.name) != sequence_hash:
                duplicate_strains.add(sequence.name)

            # If the current strain has been seen before, don't write
            # out its sequence again.
            continue

        sequence_hash_by_name[sequence.name] = sequence_hash
        yield sequence

    # Report names of duplicate strains with different sequences when requested.
    if len(duplicate_strains) > 0 and error_on_duplicates:
        raise DuplicateSequenceError(", ".join(duplicate_strains))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sequences", nargs="+", required=True, help="sequences to be sanitized")
    parser.add_argument("--strip-prefixes", nargs="+", help="prefixes to strip from strain names in the sequences")
    parser.add_argument('--error-on-duplicate-strains', action="store_true", help="exit with an error when the same strain is detected multiple times with different sequences. By default, use the first occurrence of each duplicated sequence.")
    parser.add_argument("--output", required=True, help="sanitized sequences")

    args = parser.parse_args()

    sequence_files = []
    tar_files = []
    for sequence_filename in args.sequences:
        # If the input is a tarball, try to find a sequence file inside the
        # archive.
        if ".tar" in Path(sequence_filename).suffixes:
            try:
                sequence_file, tar_file = stream_tar_file_contents(
                    sequence_filename,
                    "sequences"
                )
                sequence_files.append(sequence_file)
                tar_files.append(tar_file)
            except FileNotFoundError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                sys.exit(1)
        else:
            sequence_files.append(sequence_filename)

    # Replace whitespace and everything following pipes with nothing.
    pattern = "( )|(\|.*)"
    if args.strip_prefixes:
        prefixes = "|".join(args.strip_prefixes)
        pattern = f"^({prefixes})|{pattern}"

    with open_file(args.output, "w", threads=1) as output_handle:
        # In order to prefer the latter files, we have to reverse the order of
        # the files.
        sequences = read_sequences(*reversed(sequence_files))
        renamed_sequences = rename_sequences(sequences, pattern)
        deduplicated_sequences = drop_duplicate_sequences(
            renamed_sequences,
            args.error_on_duplicate_strains
        )

        try:
            for sequence in deduplicated_sequences:
                write_sequences(sequence, output_handle)
        except DuplicateSequenceError as error:
            print(
                f"ERROR: The following strains have duplicate sequences: {error}",
                file=sys.stderr
            )
            sys.exit(1)

    # Clean up temporary directory and files that came from a tarball.
    for tar_file in tar_files:
        tar_file.close()
