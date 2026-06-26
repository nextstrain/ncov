"""
Propose SARS-CoV-2 Pango lineages that have crossed Nextstrain clade-designation
criteria, assembling all evidence into a dossier for review.

This script only *gathers evidence*; it never edits defaults/. A human (or Claude
Code) reads results/clade_candidates.md, decides which lineage to elevate and where
to draw the demarcation, and then writes the clades.tsv / clade_display_names.yml /
clade_hierarchy.tsv / clade_emergence_dates.tsv / color_ordering.tsv edits.

Designation criteria (see https://nextstrain.org/blog/2022-04-29-SARS-CoV-2-clade-naming-2022):
a candidate Pango lineage carries >=1 spike mutation relative to its parent Nextstrain
clade and has risen to >30% regional or >20% global frequency.

Data sources (all open / GISAID-free):
  - Nextclade SARS-CoV-2 reference tree -> clean lineage-defining mutations + parent clade
  - Nextstrain open 6m builds (global/north-america/europe) tip-frequencies -> headline frequency
  - LAPIS open API (lapis.cov-spectrum.org) -> exact per-lineage regional frequency
  - forecasts-ncov open MLR -> growth advantage (fitness) + modeled/forecast frequency

Because open data is sparse and ~2 months lagged in 2026, frequencies are approximate;
the dossier surfaces every number alongside its source so they can be cross-checked.
"""

import argparse
import gzip
import json
import os
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NEXTCLADE_TREE_URL = (
    "https://nextstrain.org/charon/getDataset"
    "?prefix=/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs"
)
CHARON = "https://nextstrain.org/charon"
# Open run analyzes these core builds (dataset prefix /ncov/open/<region>/6m).
OPEN_REGIONS = ["global", "north-america", "europe"]
GISAID_GROUP_DEFAULT = "blab"
MLR_URLS = {
    "open": "https://data.nextstrain.org/files/workflows/forecasts-ncov/open"
            "/pango_lineages/global/mlr/latest_results.json",
    "gisaid": "https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid"
              "/pango_lineages/global/mlr/latest_results.json",
}
# LAPIS open instance needs no auth; the GISAID instance requires API access we no
# longer have, so LAPIS is skipped for GISAID runs.
LAPIS_URL = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"

# Region slug -> `region` node-attribute value. Used to region-filter a regional
# build's frequency (its tips are focal+contextual) and as the LAPIS region name.
REGION_NAMES = {
    "north-america": "North America",
    "south-america": "South America",
    "europe": "Europe",
    "asia": "Asia",
    "africa": "Africa",
    "oceania": "Oceania",
}

# Designation thresholds (fractions).
GLOBAL_THRESHOLD = 0.20
REGIONAL_THRESHOLD = 0.30
# Candidates below this everywhere are not worth surfacing at all.
ENUMERATION_FLOOR = 0.05
# LAPIS recent window (days) and small-sample warning levels.
LAPIS_WINDOW_DAYS = 180
LOW_DENOMINATOR = 100      # LAPIS total-count warning
LOW_REGION_TIPS = 30       # effective regional sequences backing a latest-pivot frequency

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# Fetch helpers (gzip-aware, with a simple on-disk cache)
# ---------------------------------------------------------------------------

def _http_get(url):
    req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        raw = resp.read()
    if resp.headers.get("Content-Encoding") == "gzip" or raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw


def fetch_json(url, cache_dir=None, cache_name=None):
    """GET + parse JSON, optionally caching the raw bytes under cache_dir."""
    cache_path = None
    if cache_dir and cache_name:
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = os.path.join(cache_dir, cache_name)
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as fh:
                return json.loads(fh.read())
    raw = _http_get(url)
    if cache_path:
        with open(cache_path, "wb") as fh:
            fh.write(raw)
    return json.loads(raw)


def charon_url(prefix, sidecar=None):
    """Build a nextstrain.org charon getDataset URL for a dataset path."""
    params = {"prefix": prefix}
    if sidecar:
        params["type"] = sidecar
    return f"{CHARON}/getDataset?" + urllib.parse.urlencode(params)


def _cache_key(prefix):
    return prefix.strip("/").replace("/", "_")


def fetch_build(prefix, cache_dir, refresh=False):
    """Fetch an auspice build (main tree + tip-frequencies) by dataset path -> Build."""
    key = _cache_key(prefix)
    main_cache = None if refresh else f"{key}.json"
    tf_cache = None if refresh else f"{key}_tf.json"
    tree = fetch_json(charon_url(prefix), cache_dir, main_cache)
    tf = fetch_json(charon_url(prefix, "tip-frequencies"), cache_dir, tf_cache)
    return Build(tree, tf)


def discover_gisaid_prefixes(group, pin_date=None):
    """
    Find this group's GISAID 6m datasets via charon getAvailable.

    Returns (date, {region_slug: dataset_prefix}). Picks one consistent date for all
    regions: `pin_date` if given, else the latest date the global build is available.
    Regions present at that date are returned; any missing region is the caller's to warn on.
    """
    avail = fetch_json(f"{CHARON}/getAvailable?prefix=/groups/{group}")
    # request paths look like 'groups/<group>/ncov/gisaid/<region>/6m/<date>'
    import re
    pat = re.compile(rf"^groups/{re.escape(group)}/ncov/gisaid/([a-z-]+)/6m/(\d{{4}}-\d{{2}}-\d{{2}})$")
    by_region = defaultdict(dict)  # region -> {date: prefix}
    for entry in avail.get("datasets", []):
        req = entry.get("request", "").strip("/")
        m = pat.match(req)
        if m:
            region, date = m.group(1), m.group(2)
            by_region[region][date] = "/" + req
    if not by_region:
        raise RuntimeError(f"No GISAID 6m datasets found under /groups/{group}")
    date = pin_date or max(by_region.get("global", {}), default=None)
    if date is None:
        raise RuntimeError(f"No global GISAID 6m build found under /groups/{group}")
    prefixes = {region: dates[date] for region, dates in by_region.items() if date in dates}
    return date, prefixes


def lapis_count(lineage=None, region=None, date_from=None):
    """Return the open-data sample count matching the given filters (or None on error)."""
    params = {}
    if lineage:
        params["nextcladePangoLineage"] = lineage
    if region:
        params["region"] = region
    if date_from:
        params["dateFrom"] = date_from
    url = LAPIS_URL + "?" + urllib.parse.urlencode(params)
    try:
        data = json.loads(_http_get(url))["data"]
        return sum(row["count"] for row in data)
    except Exception as exc:  # network / API hiccup -> degrade gracefully
        print(f"  ! LAPIS query failed ({exc})", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Nextclade reference tree
# ---------------------------------------------------------------------------

def _attr(node, key):
    val = node.get("node_attrs", {}).get(key)
    return val.get("value") if isinstance(val, dict) else None


class RefTree:
    """The Nextclade reference tree: clean lineage -> mutations + parent clade."""

    def __init__(self, tree_json):
        self.nodes = []
        self.parent = {}
        self.cds = {}  # gene -> (start, end, strand) for + strand single-segment CDSes
        for gene, info in tree_json.get("meta", {}).get("genome_annotations", {}).items():
            if gene == "nuc" or "segments" in info:
                continue
            start, end, strand = info.get("start"), info.get("end"), info.get("strand", "+")
            if isinstance(start, int) and isinstance(end, int):
                self.cds[gene] = (start, end, strand)
        root = tree_json["tree"]
        self._walk(root, None)
        # Origin node of each lineage = shallowest node whose pango differs from its parent's.
        self.entry = {}
        for n in self.nodes:
            lineage = _attr(n, "Nextclade_pango")
            if lineage is None:
                continue
            par = self.parent[id(n)]
            if par is None or _attr(par, "Nextclade_pango") != lineage:
                self.entry.setdefault(lineage, n)

    def _walk(self, node, par):
        self.nodes.append(node)
        self.parent[id(node)] = par
        for child in node.get("children", []):
            self._walk(child, node)

    def lineages(self):
        return set(self.entry.keys())

    def clade_of(self, lineage):
        node = self.entry.get(lineage)
        return _attr(node, "clade_nextstrain") if node is not None else None

    def who_of(self, lineage):
        node = self.entry.get(lineage)
        return _attr(node, "clade_who") if node is not None else None

    def designation_date(self, lineage):
        node = self.entry.get(lineage)
        return _attr(node, "designation_date") if node is not None else None

    def parent_lineage(self, lineage):
        node = self.entry.get(lineage)
        if node is None:
            return None
        cur = self.parent[id(node)]
        while cur is not None:
            par_lineage = _attr(cur, "Nextclade_pango")
            if par_lineage != lineage:
                return par_lineage
            cur = self.parent[id(cur)]
        return None

    def child_lineages(self, lineage):
        """Immediate descendant lineages (one Pango step below)."""
        node = self.entry.get(lineage)
        if node is None:
            return []
        out, stack = [], list(node.get("children", []))
        while stack:
            cur = stack.pop()
            cur_lineage = _attr(cur, "Nextclade_pango")
            if cur_lineage != lineage:
                if cur_lineage is not None:
                    out.append(cur_lineage)
            else:
                stack.extend(cur.get("children", []))
        # de-dup preserving order
        seen, uniq = set(), []
        for lin in out:
            if lin not in seen:
                seen.add(lin)
                uniq.append(lin)
        return uniq

    def gene_codon(self, site):
        """Map a 1-based nuc site to (gene, codon_number) for the CDS containing it."""
        for gene, (start, end, strand) in self.cds.items():
            if strand == "+" and start <= site <= end:
                return gene, (site - start) // 3 + 1
        return None, None

    def descendants(self, lineage):
        """All lineage labels in the subtree rooted at this lineage's origin (incl. self)."""
        node = self.entry.get(lineage)
        if node is None:
            return {lineage}
        labels, stack = set(), [node]
        while stack:
            cur = stack.pop()
            lab = _attr(cur, "Nextclade_pango")
            if lab:
                labels.add(lab)
            stack.extend(cur.get("children", []))
        return labels

    def _branch(self, node):
        muts = node.get("branch_attrs", {}).get("mutations", {})
        return {"lineage": _attr(node, "Nextclade_pango"),
                "nuc": list(muts.get("nuc", [])), "S": list(muts.get("S", []))}

    def defining_path(self, lineage, designated_lineages=frozenset(), lineage_to_clade=None):
        """
        Mutations distinguishing this lineage from the clade it belongs to.

        The parent clade is the nearest ancestor that is a designated clade lineage —
        preferring a locally-designated lineage (from clades.tsv/display_names) the
        reference tree may not encode yet, and otherwise falling back to the tree's
        own ``clade_nextstrain`` boundary. We collect branch mutations from the
        lineage's origin up to (but not including) the parent clade's origin; those
        define the parent clade itself and are inherited via the ``clade <parent>`` row.

        Returns (parent_clade, branches) where branches is a top-down list of
        per-step dicts {lineage, nuc, S}; their union is the net set distinguishing
        the lineage from its parent clade and a human picks the characteristic subset.
        """
        node = self.entry.get(lineage)
        if node is None:
            return None, []
        lineage_to_clade = lineage_to_clade or {}

        # Prefer the nearest ancestor lineage that is a designated clade.
        anc = None
        cur = self.parent[id(node)]
        while cur is not None:
            lab = _attr(cur, "Nextclade_pango")
            if lab != lineage and lab in designated_lineages:
                anc = lab
                break
            cur = self.parent[id(cur)]
        if anc is not None:
            stop = self.entry.get(anc)
            parent_clade = lineage_to_clade.get(anc) or _attr(stop, "clade_nextstrain")
            chain = []
            cur = node
            while cur is not None and cur is not stop:
                chain.append(self._branch(cur))
                cur = self.parent[id(cur)]
            chain.reverse()
            return parent_clade, chain

        # Fallback: the reference tree's own clade boundary.
        clade = _attr(node, "clade_nextstrain")
        chain = []
        cur = node
        while cur is not None and _attr(cur, "clade_nextstrain") == clade:
            par = self.parent[id(cur)]
            if par is None or _attr(par, "clade_nextstrain") != clade:
                break  # cur is the clade origin; its branch defines the parent clade
            chain.append(self._branch(cur))
            cur = par
        chain.reverse()
        return clade, chain


def net_nuc(branches):
    out = []
    for b in branches:
        out.extend(b["nuc"])
    return out


def net_spike(branches):
    out = []
    for b in branches:
        out.extend(b["S"])
    return out


def annotate_nuc_rows(ref, branches):
    """
    Annotate each net nuc mutation with the spike amino-acid change it causes (if any).

    Returns a list of {site, alt, nuc, spike_aa} dicts. spike_aa is e.g. 'S:E96D' when
    the nucleotide sits in spike and lines up with a tracked S substitution on the path,
    else None. The spike-linked rows are the natural characteristic subset for clades.tsv.
    """
    import re
    s_by_codon = {}
    for b in branches:
        for aa in b["S"]:
            m = re.search(r"\d+", aa)
            if m:
                s_by_codon[int(m.group())] = aa
    rows = []
    for b in branches:
        for nuc in b["nuc"]:
            site = int(nuc[1:-1])
            gene, codon = ref.gene_codon(site)
            spike_aa = f"S:{s_by_codon[codon]}" if gene == "S" and codon in s_by_codon else None
            rows.append({"site": str(site), "alt": nuc[-1], "nuc": nuc, "spike_aa": spike_aa})
    return rows


def decimal_to_date(dec):
    """Decimal year (e.g. 2025.62) -> datetime.date (approximate)."""
    import datetime
    if dec is None:
        return None
    year = int(dec)
    start = datetime.date(year, 1, 1)
    return start + datetime.timedelta(days=int(round((dec - year) * 365)))


# ---------------------------------------------------------------------------
# Repo clade state
# ---------------------------------------------------------------------------

def load_designated(clades_tsv, display_yml):
    """
    Return (clade_names, designated_lineages, lineage_to_clade) from the repo's
    default files. lineage_to_clade maps each designated Pango lineage to its clade
    name (e.g. 'XFG.1.1' -> '26A'), so candidates can be anchored to locally-added
    clades the Nextclade reference tree does not yet know about.
    """
    clade_names = set()
    with open(clades_tsv) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                clade_names.add(parts[0])

    designated_lineages = set()
    lineage_to_clade = {}
    if os.path.exists(display_yml):
        with open(display_yml) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or ":" not in line:
                    continue
                name, _, disp = line.partition(":")
                name = name.strip()
                disp = disp.strip()
                if "(" in disp and ")" in disp:
                    lineage = disp[disp.index("(") + 1:disp.index(")")].strip()
                    if lineage:
                        designated_lineages.add(lineage)
                        lineage_to_clade[lineage] = name
    return clade_names, designated_lineages, lineage_to_clade


def next_clade_name(clade_names, year):
    """First unused letter for the designation year, e.g. year=2026 -> '26A'."""
    yy = f"{year % 100:02d}"
    used = {n[2:] for n in clade_names if n[:2] == yy and len(n) == 3}
    for code in range(ord("A"), ord("Z") + 1):
        if chr(code) not in used:
            return yy + chr(code)
    return yy + "Z"  # exhausted; caller should reconsider the year


# ---------------------------------------------------------------------------
# Frequency from Nextstrain builds
# ---------------------------------------------------------------------------

class Build:
    """A Nextstrain auspice build + its tip-frequencies, for label-based frequency."""

    def __init__(self, tree_json, tipfreq_json):
        self.pivots = tipfreq_json["pivots"]
        self.freqs = {k: v["frequencies"] for k, v in tipfreq_json.items()
                      if isinstance(v, dict) and "frequencies" in v}
        self.tip_lineage = {}
        self.tip_region = {}
        self.tip_numdate = {}
        self._collect(tree_json["tree"])

    def _collect(self, node):
        if node.get("children"):
            for child in node["children"]:
                self._collect(child)
        else:
            name = node.get("name")
            self.tip_lineage[name] = _attr(node, "Nextclade_pango")
            self.tip_region[name] = _attr(node, "region")
            self.tip_numdate[name] = _attr(node, "num_date")

    def earliest_date(self, member_lineages):
        """Earliest tip collection date (decimal year) among member lineages, or None."""
        members = set(member_lineages)
        dates = [self.tip_numdate[name] for name, lin in self.tip_lineage.items()
                 if lin in members and self.tip_numdate.get(name) is not None]
        return min(dates) if dates else None

    def trajectory(self, member_lineages, region=None):
        """Summed frequency per pivot for tips in member_lineages (optionally one region)."""
        totals = [0.0] * len(self.pivots)
        denom = [0.0] * len(self.pivots)
        members = set(member_lineages)
        for name, series in self.freqs.items():
            if region is not None and self.tip_region.get(name) != region:
                continue
            for i, val in enumerate(series):
                denom[i] += val
            if self.tip_lineage.get(name) in members:
                for i, val in enumerate(series):
                    totals[i] += val
        return [t / d if d else 0.0 for t, d in zip(totals, denom)], denom

    def current_freq(self, member_lineages, region=None):
        traj, _ = self.trajectory(member_lineages, region)
        return traj[-1] if traj else 0.0

    def all_lineages(self):
        return set(filter(None, self.tip_lineage.values()))

    def region_effective_n(self, region):
        """Effective # of regional sequences backing the latest-pivot frequency.

        The KDE concentrates the final pivot on the few most-recent tips, so a plain
        count overstates support. Kish effective sample size (Σw)²/Σw² over the
        region's last-pivot tip frequencies measures how many sequences actually back
        the headline number — the signal that should drive the small-sample warning.
        """
        w = [self.freqs[name][-1] for name, r in self.tip_region.items()
             if r == region and name in self.freqs]
        s1, s2 = sum(w), sum(x * x for x in w)
        return round(s1 * s1 / s2) if s2 > 0 else 0


def trajectory_slope(traj, window=6):
    """Recent change in frequency over the last `window` pivots (model-free dynamics)."""
    if len(traj) < 2:
        return 0.0
    window = min(window, len(traj) - 1)
    return traj[-1] - traj[-1 - window]


# ---------------------------------------------------------------------------
# MLR growth advantage + modeled frequency
# ---------------------------------------------------------------------------

def load_mlr(url, cache_dir, cache_name="mlr.json"):
    """variant -> {location -> {'ga', 'freq', 'freq_forecast'}} from MLR results."""
    try:
        payload = fetch_json(url, cache_dir, cache_name)
    except Exception as exc:
        print(f"  ! MLR fetch failed ({exc})", file=sys.stderr)
        return {}
    out = defaultdict(lambda: defaultdict(dict))
    latest_freq_date = defaultdict(lambda: defaultdict(lambda: ("", None)))
    for entry in payload.get("data", []):
        if entry.get("ps") != "median":
            continue
        variant = entry.get("variant")
        loc = entry.get("location")
        site = entry.get("site")
        if variant is None or loc is None:
            continue
        if site == "ga":
            out[variant][loc]["ga"] = entry.get("value")
        elif site in ("freq", "freq_forecast"):
            date = entry.get("date", "")
            prev_date, _ = latest_freq_date[variant][(loc, site)]
            if date >= prev_date:  # keep the most recent modeled point
                latest_freq_date[variant][(loc, site)] = (date, entry.get("value"))
    for variant, locsite in latest_freq_date.items():
        for (loc, site), (_, value) in locsite.items():
            out[variant][loc][site] = value
    return {v: dict(d) for v, d in out.items()}


def hier_ga(mlr, lineage):
    return mlr.get(lineage, {}).get("hierarchical", {}).get("ga")


# ---------------------------------------------------------------------------
# Candidate detection
# ---------------------------------------------------------------------------

def build_candidate(lineage, ref, builds, mlr, designated_lineages, lineage_to_clade,
                    use_lapis, lapis_min, date_from):
    """Assemble the full evidence record for one candidate lineage, or None to skip."""
    descendants = ref.descendants(lineage)
    # Skip ancestral catch-alls: a lineage that already contains a designated clade
    # is an ancestor you'd refine, not a new clade in its own right.
    if any(d in designated_lineages for d in descendants if d != lineage):
        return None

    parent_clade, branches = ref.defining_path(lineage, designated_lineages, lineage_to_clade)
    if parent_clade is None or not branches:
        return None
    spike = net_spike(branches)
    nuc = net_nuc(branches)
    nuc_rows = annotate_nuc_rows(ref, branches)
    is_recombinant = lineage.startswith("X")

    members = set(descendants)
    for build in builds.values():  # include build labels newer than the ref tree
        members |= {l for l in build.all_lineages()
                    if l == lineage or l.startswith(lineage + ".")}

    freq = {"build": {}, "lapis": {}, "mlr": {}}
    gb = builds["global"]
    traj_global, _ = gb.trajectory(members)  # global = no region filter
    freq["build"]["global"] = traj_global[-1] if traj_global else 0.0
    # Each regional build is focal+contextual, so a region's frequency is taken over only
    # the tips whose `region` attribute matches (the ?f_region=<Region> equivalent).
    regional_keys = [k for k in builds if k != "global"]
    for rk in regional_keys:
        freq["build"][rk] = builds[rk].current_freq(members, region=REGION_NAMES.get(rk))

    g = freq["build"]["global"]
    regional_vals = [freq["build"][rk] for rk in regional_keys]
    peak_build = max([g] + regional_vals)
    if peak_build < ENUMERATION_FLOOR:
        return None
    passes_build = g > GLOBAL_THRESHOLD or any(v > REGIONAL_THRESHOLD for v in regional_vals)

    # LAPIS is the exact read but costs queries — only run it for the shortlist.
    if use_lapis and (passes_build or peak_build >= lapis_min):
        for rk in regional_keys:
            region_val = REGION_NAMES.get(rk)
            if region_val is None:
                continue
            lin_count = lapis_count(lineage + "*", region_val, date_from)
            total = lapis_count(None, region_val, date_from)
            if lin_count is not None and total:
                freq["lapis"][rk] = {
                    "freq": lin_count / total, "n": lin_count, "total": total}
        g_lin = lapis_count(lineage + "*", None, date_from)
        g_total = lapis_count(None, None, date_from)
        if g_lin is not None and g_total:
            freq["lapis"]["global"] = {"freq": g_lin / g_total, "n": g_lin, "total": g_total}

    for loc, vals in mlr.get(lineage, {}).items():
        freq["mlr"][loc] = vals

    earliest_dec = gb.earliest_date(members)
    earliest_date = decimal_to_date(earliest_dec)
    first_sequence = f"{earliest_date.year:04d}-{earliest_date.month:02d}-01" if earliest_date else None

    lapis_regional = [v["freq"] for k, v in freq["lapis"].items() if k != "global"]
    lapis_g = freq["lapis"].get("global", {}).get("freq", 0.0)
    passes_lapis = (lapis_g > GLOBAL_THRESHOLD
                    or any(v > REGIONAL_THRESHOLD for v in lapis_regional))
    passes = passes_build or passes_lapis
    peak_freq = max([peak_build, lapis_g] + lapis_regional)

    return {
        "lineage": lineage,
        "parent_clade": parent_clade,
        "who": ref.who_of(lineage),
        "designation_date": ref.designation_date(lineage),
        "earliest_date": earliest_dec,
        "first_sequence": first_sequence,
        "is_recombinant": is_recombinant,
        "spike_mutations": spike,
        "net_nuc": nuc,
        "nuc_rows": nuc_rows,
        "spike_nuc_rows": [r for r in nuc_rows if r["spike_aa"]],
        "branches": branches,
        "frequency": freq,
        "trajectory_global": traj_global,
        "slope_global": trajectory_slope(traj_global),
        "mlr_ga": hier_ga(mlr, lineage),
        "peak_freq": peak_freq,
        "has_spike": len(spike) > 0,
        "passes_threshold": passes,
        "passes_build": passes_build,
        "flagged": passes and len(spike) > 0,
    }


def demarcation_chain(lineage, ref, builds, mlr):
    """Parent lineage + candidate + child lineages, each with freq + Δga, for the split call."""
    chain = []
    options = []
    parent = ref.parent_lineage(lineage)
    if parent:
        options.append(("parent", parent))
    options.append(("candidate", lineage))
    for child in ref.child_lineages(lineage)[:6]:
        options.append(("child", child))

    global_build = builds.get("global")
    for role, lin in options:
        members = ref.descendants(lin)
        cur = global_build.current_freq(members) if global_build else None
        traj, _ = global_build.trajectory(members) if global_build else ([], [])
        ga = hier_ga(mlr, lin)
        par = ref.parent_lineage(lin)
        delta_ga = None
        if ga is not None and par is not None and hier_ga(mlr, par) is not None:
            delta_ga = ga - hier_ga(mlr, par)
        chain.append({
            "role": role,
            "lineage": lin,
            "parent_lineage": par,
            "global_freq": cur,
            "slope_global": trajectory_slope(traj) if traj else None,
            "mlr_ga": ga,
            "delta_ga": delta_ga,
        })
    return chain


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def pct(x):
    return "n/a" if x is None else f"{100 * x:.0f}%"


def _region_freq_line(cand, meta):
    """global X | <Region> Y (⚠low-n) | ... for whatever regions were analyzed."""
    b = cand["frequency"]["build"]
    region_tips = meta.get("region_tips", {})
    parts = [f"global {pct(b.get('global'))}"]
    for rk in meta.get("regions", []):
        n = region_tips.get(rk, 0)
        warn = f" ⚠low-n={n}" if n < LOW_REGION_TIPS else ""
        parts.append(f"{REGION_NAMES.get(rk, rk)} {pct(b.get(rk))}{warn}")
    return " | ".join(parts)


def write_markdown(path, flagged, watch, suggested_names, meta):
    source = meta.get("source", "open")
    date = meta.get("gisaid_date")
    src = f"**{source}**" + (f" ({date})" if date else "")
    region_list = ", ".join(["global"] + [REGION_NAMES.get(r, r) for r in meta.get("regions", [])])
    lines = []
    lines.append("# Clade designation candidates\n")
    lines.append(f"_Generated from {src} data. Ref tree updated {meta.get('ref_updated','?')}; "
                 f"global build latest pivot {meta.get('latest_pivot','?')}._\n")
    lines.append(f"Regions analyzed: {region_list}.\n")
    lines.append("Criteria: >=1 spike mutation vs parent clade AND (>20% global OR >30% in any "
                 "region). Regional frequency filters tips to that region; data is sparse/lagged, so "
                 "treat percentages as approximate and cross-check the live resources.\n")
    lines.append("Suggested names use the **designation year** (the current year), not the "
                 "lineage's emergence year — so the next name is the first unused letter for "
                 "this year.\n")

    if not flagged:
        lines.append("\n**No lineages currently pass the designation threshold.**\n")
    for cand in flagged:
        _render_candidate(lines, cand, suggested_names.get(cand["lineage"]), meta, header="## ")

    lines.append("\n---\n## Watch list (rising, not yet over threshold)\n")
    if not watch:
        lines.append("_None above the enumeration floor._\n")
    for cand in watch:
        lines.append(
            f"- **{cand['lineage']}** (parent clade {cand['parent_clade']}, "
            f"{len(cand['spike_mutations'])} spike mut) — {_region_freq_line(cand, meta)}; "
            f"MLR ga {cand['mlr_ga'] if cand['mlr_ga'] is not None else 'n/a'}"
        )

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _render_candidate(lines, cand, suggested, meta, header="## "):
    lin = cand["lineage"]
    lines.append(f"\n{header}{lin}  →  suggested name **{suggested}** "
                 f"(parent clade {cand['parent_clade']})")
    if cand.get("alternatives"):
        lines.append(f"\n> Alternative demarcation of the same sweep — choose ONE of "
                     f"{lin} or {', '.join(cand['alternatives'])} to carry {suggested}.")
    if cand["is_recombinant"]:
        lines.append("\n> ⚠ Recombinant (X*) lineage — sanity-check the breakpoint and that the "
                     "defining mutations are inherited cleanly, not split across parents.")
    lines.append("\n**Frequency (Nextstrain 6m build tip-frequencies — headline):**")
    lines.append(f"- {_region_freq_line(cand, meta)}")
    if cand["frequency"]["lapis"]:
        parts = []
        for rk in ["global"] + meta.get("regions", []):
            if rk in cand["frequency"]["lapis"]:
                v = cand["frequency"]["lapis"][rk]
                warn = "  ⚠low-n" if v["total"] < LOW_DENOMINATOR else ""
                label = "global" if rk == "global" else REGION_NAMES.get(rk, rk)
                parts.append(f"{label} {pct(v['freq'])} ({v['n']}/{v['total']}{warn})")
        lines.append("\n**Frequency (LAPIS open, exact, 180d):** " + " | ".join(parts))
    if cand["frequency"]["mlr"]:
        parts = []
        for loc, vals in cand["frequency"]["mlr"].items():
            seg = []
            if vals.get("ga") is not None:
                seg.append(f"ga={vals['ga']}")
            if vals.get("freq") is not None:
                seg.append(f"freq={pct(vals['freq'])}")
            if vals.get("freq_forecast") is not None:
                seg.append(f"forecast={pct(vals['freq_forecast'])}")
            if seg:
                parts.append(f"{loc}: " + ", ".join(seg))
        if parts:
            lines.append("\n**MLR fitness/frequency:** " + " | ".join(parts))
    lines.append(f"\n**Spike mutations vs parent clade ({len(cand['spike_mutations'])}):** "
                 + (", ".join(cand["spike_mutations"]) or "none"))
    lines.append(f"\n**Global frequency trend (last 6 pivots):** "
                 f"{cand['slope_global']:+.2f}")
    lines.append(f"\n**First seen:** {cand.get('first_sequence') or 'unknown'}"
                 f"  ·  WHO: {cand.get('who') or 'n/a'}")

    chain = cand.get("demarcation", [])
    if chain:
        lines.append("\n**Where to draw the line (prefer larger fitness differential):**")
        lines.append("\n| role | lineage | global freq | trend | MLR ga | Δga vs parent |")
        lines.append("|---|---|---|---|---|---|")
        for c in chain:
            lines.append(
                f"| {c['role']} | {c['lineage']} | {pct(c['global_freq'])} | "
                f"{('%+.2f' % c['slope_global']) if c['slope_global'] is not None else 'n/a'} | "
                f"{c['mlr_ga'] if c['mlr_ga'] is not None else 'n/a'} | "
                f"{('%+.2f' % c['delta_ga']) if c['delta_ga'] is not None else 'n/a'} |"
            )

    lineage = cand["lineage"]
    who = cand.get("who") or "Omicron"
    emergence = cand.get("first_sequence") or "YYYY-MM-01"

    def fmt_rows(rows):
        return "\n".join(
            f"{suggested}\tnuc\t{r['site']}\t{r['alt']}"
            + (f"\t# {r['spike_aa']}" if r["spike_aa"] else "")
            for r in rows)

    lines.append("\n**Proposed file edits** (the clades.tsv block lists ALL net mutations vs "
                 "the parent clade; trim to a characteristic handful — the spike-linked rows, "
                 "flagged below, are the natural choice):")
    lines.append("```")
    lines.append(f"# defaults/clades.tsv  (net mutations vs parent clade {cand['parent_clade']})")
    lines.append(f"{suggested}\tclade\t{cand['parent_clade']}")
    lines.append(fmt_rows(cand["nuc_rows"]))
    if cand["spike_nuc_rows"]:
        lines.append(f"#   suggested characteristic subset (spike-linked): "
                     + ", ".join(r["spike_aa"] for r in cand["spike_nuc_rows"]))
    lines.append(f"\n# defaults/clade_display_names.yml")
    lines.append(f"{suggested}: {suggested} ({lineage})")
    lines.append(f"\n# defaults/clade_hierarchy.tsv")
    lines.append(f"{suggested}\t{cand['parent_clade']}\t{who}")
    lines.append(f"\n# defaults/clade_emergence_dates.tsv")
    lines.append(f"{suggested}\t{emergence}")
    lines.append(f"\n# defaults/color_ordering.tsv")
    lines.append(f"clade_membership\t{suggested} ({lineage})")
    lines.append("```")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-source", choices=["open", "gisaid"], default="open",
                        help="open: fetch live /ncov/open/{global,north-america,europe}/6m. "
                             "gisaid: auto-discover the group's dated /groups/<group>/ncov/gisaid/<region>/6m builds.")
    parser.add_argument("--group", default=GISAID_GROUP_DEFAULT, help="Nextstrain group for GISAID builds")
    parser.add_argument("--gisaid-date", default=None, help="Pin a GISAID upload date (default: latest available)")
    parser.add_argument("--gisaid-url-template", default=None,
                        help="Override discovery with an explicit dataset-path template containing '{region}', "
                             "e.g. /groups/blab/ncov/gisaid/{region}/6m/2026-06-25")
    parser.add_argument("--clades", default=os.path.join(REPO_ROOT, "defaults", "clades.tsv"))
    parser.add_argument("--display-names", default=os.path.join(REPO_ROOT, "defaults", "clade_display_names.yml"))
    parser.add_argument("--cache-dir", default=os.path.join(REPO_ROOT, "results", "clade_cache"))
    parser.add_argument("--output-md", default=os.path.join(REPO_ROOT, "results", "clade_candidates.md"))
    parser.add_argument("--output-json", default=os.path.join(REPO_ROOT, "results", "clade_candidates.json"))
    parser.add_argument("--date-from", default=None, help="LAPIS window start (default: 180 days before latest pivot)")
    parser.add_argument("--designation-year", type=int, default=None, help="Year for clade names (default: current year, i.e. the designation year)")
    parser.add_argument("--lapis-min", type=float, default=0.15, help="Only run LAPIS for candidates at/above this build frequency")
    parser.add_argument("--skip-lapis", action="store_true", help="Skip exact LAPIS frequency queries")
    parser.add_argument("--skip-regional", action="store_true", help="Analyze the global build only")
    parser.add_argument("--skip-mlr", action="store_true", help="Skip forecasts-ncov MLR fitness")
    parser.add_argument("--refresh", action="store_true", help="Ignore cached downloads")
    args = parser.parse_args()

    cache_dir = args.cache_dir
    if args.refresh and os.path.isdir(cache_dir):
        for f in os.listdir(cache_dir):
            os.remove(os.path.join(cache_dir, f))

    print("Loading Nextclade reference tree ...", file=sys.stderr)
    ref_json = fetch_json(NEXTCLADE_TREE_URL, cache_dir, "nextclade_ref.json")
    ref = RefTree(ref_json)
    ref_updated = ref_json.get("meta", {}).get("updated", "?")

    source = args.data_source
    gisaid_date = None
    if source == "open":
        prefixes = {r: f"/ncov/open/{r}/6m" for r in OPEN_REGIONS}
    elif args.gisaid_url_template:
        prefixes = {r: args.gisaid_url_template.format(region=r)
                    for r in ["global"] + list(REGION_NAMES)}
    else:
        print(f"Discovering GISAID builds under /groups/{args.group} ...", file=sys.stderr)
        gisaid_date, prefixes = discover_gisaid_prefixes(args.group, args.gisaid_date)
        print(f"  date {gisaid_date}: {', '.join(sorted(prefixes))}", file=sys.stderr)

    if "global" not in prefixes:
        sys.exit("ERROR: no global build available for this data source.")
    if args.skip_regional:
        prefixes = {"global": prefixes["global"]}

    builds = {}
    for region, prefix in prefixes.items():
        try:
            print(f"Fetching {source} {region} build ...", file=sys.stderr)
            builds[region] = fetch_build(prefix, cache_dir)
        except Exception as exc:
            if region == "global":
                sys.exit(f"ERROR: could not load global build {prefix}: {exc}")
            print(f"  ! skipping {region} build ({exc})", file=sys.stderr)
    latest_pivot = builds["global"].pivots[-1]

    mlr = {} if args.skip_mlr else load_mlr(MLR_URLS[source], cache_dir, f"mlr_{source}.json")
    use_lapis = (source == "open") and not args.skip_lapis

    clade_names, designated_lineages, lineage_to_clade = load_designated(args.clades, args.display_names)

    date_from = args.date_from
    if date_from is None:
        # ~180 days before the latest pivot (decimal year -> ISO date, approximate).
        import datetime
        year = int(latest_pivot)
        rem = latest_pivot - year
        day_of_year = int(round(rem * 365))
        latest_date = datetime.date(year, 1, 1) + datetime.timedelta(days=day_of_year)
        date_from = (latest_date - datetime.timedelta(days=LAPIS_WINDOW_DAYS)).isoformat()

    # Enumerate observed lineages (build + MLR-tracked), drop already-designated ones.
    observed = set()
    for build in builds.values():
        observed |= build.all_lineages()
    observed |= set(mlr.keys())
    observed.discard("other")
    candidates_in = sorted(l for l in observed
                           if l in ref.lineages() and l not in designated_lineages)

    print(f"Evaluating {len(candidates_in)} non-designated lineages ...", file=sys.stderr)
    records = []
    for lineage in candidates_in:
        rec = build_candidate(lineage, ref, builds, mlr, designated_lineages, lineage_to_clade,
                              use_lapis=use_lapis, lapis_min=args.lapis_min,
                              date_from=date_from)
        if rec is not None:
            records.append(rec)

    flagged = [r for r in records if r["flagged"]]
    watch = [r for r in records if not r["flagged"]]
    flagged.sort(key=lambda r: r["peak_freq"], reverse=True)
    watch.sort(key=lambda r: r["peak_freq"], reverse=True)

    # Group flagged candidates that are nested in the tree — these are alternative
    # demarcations of the same sweeping entity, so they share one suggested name (you
    # designate exactly one of them).
    clusters = []
    for rec in flagged:  # peak-descending, so the broader cut anchors each cluster
        for cl in clusters:
            anchor = cl[0]["lineage"]
            if (rec["lineage"] in ref.descendants(anchor)
                    or anchor in ref.descendants(rec["lineage"])):
                cl.append(rec)
                break
        else:
            clusters.append([rec])

    # Suggested names + demarcation chains for flagged candidates. The year in a
    # Nextstrain clade name is the DESIGNATION year (when it is named), not the
    # lineage's emergence year — so all names use the current year.
    import datetime
    designation_year = args.designation_year or datetime.date.today().year
    suggested_names = {}
    working_names = set(clade_names)
    for cl in clusters:
        name = next_clade_name(working_names, designation_year)
        working_names.add(name)
        siblings = [r["lineage"] for r in cl]
        for rec in cl:
            suggested_names[rec["lineage"]] = name
            rec["suggested_name"] = name
            rec["alternatives"] = [l for l in siblings if l != rec["lineage"]]
            rec["demarcation"] = demarcation_chain(rec["lineage"], ref, builds, mlr)

    regional_keys = [r for r in builds if r != "global"]
    region_tips = {r: builds[r].region_effective_n(REGION_NAMES.get(r)) for r in regional_keys}
    meta = {
        "source": source,
        "gisaid_date": gisaid_date,
        "regions": regional_keys,
        "region_tips": region_tips,
        "ref_updated": ref_updated,
        "latest_pivot": round(latest_pivot, 3),
        "thresholds": {"global": GLOBAL_THRESHOLD, "regional": REGIONAL_THRESHOLD},
    }

    os.makedirs(os.path.dirname(args.output_md), exist_ok=True)
    write_markdown(args.output_md, flagged, watch, suggested_names, meta)
    with open(args.output_json, "w") as fh:
        json.dump({"meta": {**meta, "date_from": date_from},
                   "flagged": flagged, "watch": watch}, fh, indent=2)

    print(f"\nFlagged {len(flagged)} candidate(s); {len(watch)} on the watch list.", file=sys.stderr)
    print(f"Wrote {args.output_md}", file=sys.stderr)
    print(f"Wrote {args.output_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
