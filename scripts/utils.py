import lzma
from pathlib import Path
import tarfile


EXTENSION_BY_FILETYPE = {
    "metadata": ".tsv",
    "sequences": ".fasta",
}


def extract_tar_file_contents(filename, filetype):
    """Try to extract the contents of a given file type (e.g., metadata or
    sequences) from the given tar filename.

    """
    extension = EXTENSION_BY_FILETYPE[filetype]

    tar = tarfile.open(filename)
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
            elif extension == ".fasta":
                # For sequence data, handle decoding of the binary stream prior
                # to passing the data back to the caller.
                internal_file = (line.decode("utf-8") for line in internal_file)

            break

    if internal_file is None:
        tar.close()
        raise FileNotFoundError(f"Could not find a {filetype} file in '{filename}'")

    return internal_file, tar
