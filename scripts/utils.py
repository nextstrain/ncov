from pathlib import Path
import sys
import tarfile
import tempfile


EXTENSION_BY_FILETYPE = {
    "metadata": ".tsv",
    "sequences": ".fasta",
}


def extract_tar_file_contents(filename, filetype):
    """Try to extract the contents of a given file type (e.g., metadata or
    sequences) from the given tar filename.

    Parameters
    ----------
    filename : str or Path-like
        Path to the tar archive to search for the given file type.
    filetype : str
        Type of file to search for in the given tar archive based on the
        associated file extension.

    Returns
    -------
    tempfile.TemporaryDir :
        Temporary directory containing the file extracted from the tar archive.
    pathlib.Path :
        Path to the file extracted from the archive with the same name as the
        file in the original archive.

    Raises
    ------
    FileNotFoundError :
        When a file with the the requested file type's extension could not be
        found in the given tar archive.

    """
    extension = EXTENSION_BY_FILETYPE[filetype]

    with tarfile.open(filename) as tar:
        internal_member = None
        for member in tar.getmembers():
            suffixes = Path(member.name).suffixes

            if extension in suffixes:
                # By default, return the binary stream for the member file.
                internal_member = member
                break

        if internal_member is None:
            raise FileNotFoundError(f"Could not find a {filetype} file in '{filename}'")

        # Extract the internal file with its original name in the tar archive to
        # a temporary directory. This approach allows downstream processes to
        # re-read the file in multiple passes instead of making a single pass
        # through a stream.
        temporary_dir = tempfile.TemporaryDirectory()
        tar.extractall(
            temporary_dir.name,
            members=(internal_member,)
        )

        extracted_file_path = Path(temporary_dir.name) / Path(internal_member.name)
        print(f"Extracted {filetype} file from {filename} to {extracted_file_path}", file=sys.stderr)

    # Return temporary directory with the path to the extract file to allow the
    # caller to clean up this directory and to maintain a reference to this
    # directory until it is no longer needed. Python will automatically clean up
    # the temporary directory when its object is destroyed. For more details, see
    # https://docs.python.org/3/library/tempfile.html#tempfile.TemporaryDirectory
    return temporary_dir, extracted_file_path
