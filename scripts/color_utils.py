"""
Shared helpers for assigning Nextstrain colors, used by both scripts/assign-clade-colors.py
(clades + geographic traits) and scripts/assign-lineage-colors.py (Pango lineages).

Keeping the recency logic here ensures the set of "colored clades" is computed
identically in both, so clade_membership and Nextclade_pango colorings correspond.
"""

import datetime

import isodate


def relative_date(duration):
    """ISO 8601 backwards-looking duration (optional 'P' prefix, e.g. '1W'/'P1W') -> date."""
    if not duration.startswith("P"):
        duration = "P" + duration
    return datetime.date.today() - isodate.parse_duration(duration)


def date_within_window(date_str, cutoff_date):
    """True if a 'YYYY-MM-DD' date is on/after cutoff_date (uncertain 'XX' dates excluded)."""
    if "XX" in date_str:
        return False
    try:
        return datetime.datetime.strptime(date_str, "%Y-%m-%d").date() >= cutoff_date
    except ValueError:
        return False


def present_clades(clade_nodes, metadata_df, clade_recency):
    """
    Set of clade_membership values that should be colored: those with >=1 sequence
    collected within `clade_recency` (a cutoff date from relative_date), or — when
    no recency / metadata is given — every clade present in the tree. Mirrors the
    logic in scripts/assign-clade-colors.py.
    """
    if clade_recency is not None and metadata_df is not None:
        dates = dict(zip(metadata_df["strain"], metadata_df["date"]))
        return {info["clade_membership"]
                for strain, info in clade_nodes.items()
                if info.get("clade_membership") and strain in dates
                and date_within_window(str(dates[strain]), clade_recency)}
    return {info["clade_membership"]
            for info in clade_nodes.values() if info.get("clade_membership")}


def load_color_schemes(path):
    """schemes[N] = list of N colors (line N of the color-schemes TSV)."""
    schemes = {}
    with open(path) as fh:
        for i, line in enumerate(fh, start=1):
            schemes[i] = line.strip().split("\t")
    return schemes


def assign_scheme(values, schemes):
    """Colors for `values`: the len(values)-color scheme, reusing colors if too few."""
    n = len(values)
    if n == 0:
        return []
    if len(schemes) < n:  # not enough distinct colors — repeat the largest scheme
        out, remain = [], n
        while remain > 0:
            take = len(schemes) if remain > len(schemes) else remain
            out += schemes[take]
            remain -= take
        return out
    return schemes[n]
