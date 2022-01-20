"""
Shared functions to transparently handle local and remote inputs/outputs.

Snakemake itself tries to offer a similar sort of scheme-level transparency for
remote files through its ``AUTO`` provider, but it's broken.

The ``AUTO`` provider enumerates and loads all the various
``snakemake.remote.*`` modules on each call to ``AUTO.remote()``.  This slows
things down unnecessarily, though isn't a deal breaker.  However, during
enumeration, it protects against unsupported modules due to missing deps, but
*not* supported modules which fail to initialize because of missing required
credentials/config from the environment or missing required parameters to the
remote provider constructor or ``AUTO.remote()``.  This alone makes it
completely unusable.  Even with a usable ``AUTO``, we'd still need some
additional handling for local paths.

We could make a better functioning ``AUTO`` equivalent ourselves by similarly
building the scheme → remote function mapping dynamically, but it didn't seem
worth it to support use cases that aren't concrete right now (http(s), s3, gs
are all known use cases for us).  The implementation here could certainly be
further developed in that direction in the future, however.
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
    "gs": _remote_producer("GS", lambda GS: GS.RemoteProvider().remote),
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

    url = SCHEME_REMOTES[url.scheme]()(
        schemeless_url,
        stay_on_remote = stay_on_remote,
        keep_local     = keep_local,
    )

    # Remote files may be returned as lists with a single item, so we need to
    # flatten to a scalar.
    if isinstance(url, list):
        url = url[0]

    return url
