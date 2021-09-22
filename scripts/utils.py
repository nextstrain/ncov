import lzma
from pathlib import Path
import tarfile
from tempfile import NamedTemporaryFile


EXTENSION_BY_FILETYPE = {
    "metadata": ".tsv",
    "sequences": ".fasta",
}


def extract_tar_file_contents(filename, filetype):
    """Try to extract the contents of a given file type (e.g., metadata or
    sequences) from the given tar filename.

    """
    extension = EXTENSION_BY_FILETYPE[filetype]

    with tarfile.open(filename) as tar:
        internal_file = None
        for member in tar.getmembers():
            suffixes = Path(member.name).suffixes

            if extension in suffixes:
                # By default, return the binary stream for the member file.
                internal_file = tar.extractfile(member.name)

                if ".xz" in suffixes:
                    # Check for LZMA-compressed data and open these with the
                    # corresponding library.
                    internal_file = lzma.open(internal_file, "rt")

                break

        if internal_file is None:
            raise FileNotFoundError(f"Could not find a {filetype} file in '{filename}'")

        # Write the internal file out to an uncompressed temporary file. This
        # approach allows downstream processes to re-read the file in multiple
        # passes instead of making a single pass through a stream.
        with NamedTemporaryFile(delete=False) as temporary_file:
            temporary_file_path = temporary_file.name
            for line in internal_file:
                if not isinstance(line, bytes):
                    line = line.encode("utf-8")

                temporary_file.write(line)

    return temporary_file_path
