"""Shared functions to transparently handle local and remote inputs/outputs.
"""
from functools import partial
from snakemake.remote import HTTP, S3
from urllib.parse import urlsplit


SCHEME_REMOTES = {
    "https": HTTP.RemoteProvider().remote,
    "http": partial(HTTP.RemoteProvider().remote, insecure = True),
    "s3": S3.RemoteProvider().remote,
    "": lambda local_path, **kwargs: local_path,
}


def path_or_url(path_or_url, stay_on_remote = False, keep_local = False):
    """
    Wrap rule inputs and outputs which may be local paths or remote URLs.

    Automatically maps well-known URL schemes to the appropriate remote
    providers in the ``snakemake.remote.*`` classes.  The returned objects
    should be interchangeable for most rules.  For example:

    ::

        rule align:
            input: path_or_url(config["sequences"])
            output: "results/aligned.fasta"
            shell:
                '''
                augur align --sequences {input} --output {output} …
                '''

    If ``config["sequences"]`` is a URL, then Snakemake will first temporarily
    download the referenced file to a local path and then substitute that local
    path as the new value of ``input`` when running the ``shell`` block.

    See `https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html`__
    for more information on Snakemake's remote files support, including the
    supported keyword args ``keep_local`` and ``stay_on_remote``.
    """
    url = urlsplit(path_or_url)

    if url.scheme not in SCHEME_REMOTES:
        raise ValueError(f"Unsupported URL scheme: {path_or_url!r} (supported schemes are {set(SCHEME_REMOTES)})")

    return SCHEME_REMOTES[url.scheme](
        url.netloc + url.path,
        stay_on_remote = stay_on_remote,
        keep_local     = keep_local,
    )
