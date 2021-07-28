"""Shared functions to transparently handle local and remote inputs/outputs.
"""
from functools import partial
from importlib import import_module
from urllib.parse import urlsplit, urljoin


def _remote_producer(remote_name, init):
    """
    Produce a remote input function by loading the :class:`snakemake.remote`
    module named by *remote_name* and passing the module to an *init* function.

    *init* will be called with the loaded module and is expected to return a
    remote input function suitable for use inside :func:`path_or_url`.

    Returns a "producer" function which when called with no arguments returns
    the remote input function.

    The call to *init* (that is, the production of the remote input function)
    is deferred until the first time the producer is called.

    If the remote module is not available because required dependencies aren't
    installed, then a :exc:`WorkflowError` exception object is returned on the
    first call of the producer.
    """
    # Name of first argument to *init*.
    remote = None

    def producer():
        nonlocal remote
        if not remote:
            remote = init(import_module(f"snakemake.remote.{remote_name}"))
        return remote

    return producer


SCHEME_REMOTES = {
    "https": _remote_producer("HTTP", lambda HTTP: HTTP.RemoteProvider().remote),
    "http": _remote_producer("HTTP", lambda HTTP: partial(HTTP.RemoteProvider().remote, insecure = True)),
    "s3": _remote_producer("S3", lambda S3: S3.RemoteProvider().remote),
    "": lambda: lambda local_path, **kwargs: local_path,
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

    # Snakemake's providers want scheme-less URLs without even the leading
    # hierarchy root prefix (//).  urljoin() will discard the query and
    # fragment (including their leading delimiters) if they're empty.
    schemeless_url = urljoin(url.netloc + url.path, f"?{url.query}#{url.fragment}")

    return SCHEME_REMOTES[url.scheme]()(
        schemeless_url,
        stay_on_remote = stay_on_remote,
        keep_local     = keep_local,
    )
