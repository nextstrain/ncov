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
BUILD_URL = "https://data.nextstrain.org/ncov_open_{region}_6m.json"
TIP_FREQ_URL = "https://data.nextstrain.org/ncov_open_{region}_6m_tip-frequencies.json"
MLR_URL = (
    "https://data.nextstrain.org/files/workflows/forecasts-ncov/open"
    "/pango_lineages/global/mlr/latest_results.json"
)
LAPIS_URL = "https://lapis.cov-spectrum.org/open/v2/sample/aggregated"

# Region keys: build-name suffix -> LAPIS region value
REGIONS = {
    "north-america": "North America",
    "europe": "Europe",
}

# Designation thresholds (fractions).
GLOBAL_THRESHOLD = 0.20
REGIONAL_THRESHOLD = 0.30
# Candidates below this everywhere are not worth surfacing at all.
ENUMERATION_FLOOR = 0.05
# LAPIS recent window (days) and small-sample warning level.
LAPIS_WINDOW_DAYS = 180
LOW_DENOMINATOR = 100

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


def trajectory_slope(traj, window=6):
    """Recent change in frequency over the last `window` pivots (model-free dynamics)."""
    if len(traj) < 2:
        return 0.0
    window = min(window, len(traj) - 1)
    return traj[-1] - traj[-1 - window]


# ---------------------------------------------------------------------------
# MLR growth advantage + modeled frequency
# ---------------------------------------------------------------------------

def load_mlr(url, cache_dir):
    """variant -> {location -> {'ga', 'freq', 'freq_forecast'}} from open MLR results."""
    try:
        payload = fetch_json(url, cache_dir, "mlr_open.json")
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
    traj_global = []
    for region_key, build in builds.items():
        if region_key == "global":
            traj_global, _ = build.trajectory(members)
            freq["build"]["global"] = build.current_freq(members)
        else:
            freq["build"][region_key] = build.current_freq(members)

    g = freq["build"].get("global", 0.0)
    na = freq["build"].get("north-america", 0.0)
    eu = freq["build"].get("europe", 0.0)
    peak_build = max(g, na, eu)
    if peak_build < ENUMERATION_FLOOR:
        return None
    passes_build = g > GLOBAL_THRESHOLD or na > REGIONAL_THRESHOLD or eu > REGIONAL_THRESHOLD

    # LAPIS is the exact read but costs queries — only run it for the shortlist.
    if use_lapis and (passes_build or peak_build >= lapis_min):
        for region_key, region_val in REGIONS.items():
            lin_count = lapis_count(lineage + "*", region_val, date_from)
            total = lapis_count(None, region_val, date_from)
            if lin_count is not None and total:
                freq["lapis"][region_key] = {
                    "freq": lin_count / total, "n": lin_count, "total": total}
        g_lin = lapis_count(lineage + "*", None, date_from)
        g_total = lapis_count(None, None, date_from)
        if g_lin is not None and g_total:
            freq["lapis"]["global"] = {"freq": g_lin / g_total, "n": g_lin, "total": g_total}

    for loc, vals in mlr.get(lineage, {}).items():
        freq["mlr"][loc] = vals

    earliest_dec = builds["global"].earliest_date(members) if "global" in builds else None
    earliest_date = decimal_to_date(earliest_dec)
    first_sequence = f"{earliest_date.year:04d}-{earliest_date.month:02d}-01" if earliest_date else None

    lapis_na = freq["lapis"].get("north-america", {}).get("freq", 0.0)
    lapis_eu = freq["lapis"].get("europe", {}).get("freq", 0.0)
    lapis_g = freq["lapis"].get("global", {}).get("freq", 0.0)
    passes_lapis = (lapis_g > GLOBAL_THRESHOLD or lapis_na > REGIONAL_THRESHOLD
                    or lapis_eu > REGIONAL_THRESHOLD)
    passes = passes_build or passes_lapis

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
        "peak_freq": max(peak_build, lapis_na, lapis_eu, lapis_g),
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


def write_markdown(path, flagged, watch, suggested_names, meta):
    lines = []
    lines.append("# Clade designation candidates\n")
    lines.append(f"_Generated from open data. Ref tree updated {meta.get('ref_updated','?')}; "
                 f"global build latest pivot {meta.get('latest_pivot','?')}._\n")
    lines.append("Criteria: >=1 spike mutation vs parent clade AND (>20% global OR >30% "
                 "regional). Open data is sparse/lagged — treat percentages as approximate "
                 "and cross-check against the live resources.\n")
    lines.append("Suggested names use the **designation year** (the current year), not the "
                 "lineage's emergence year — so the next name is the first unused letter for "
                 "this year.\n")

    if not flagged:
        lines.append("\n**No lineages currently pass the designation threshold.**\n")
    for cand in flagged:
        _render_candidate(lines, cand, suggested_names.get(cand["lineage"]), header="## ")

    lines.append("\n---\n## Watch list (rising, not yet over threshold)\n")
    if not watch:
        lines.append("_None above the enumeration floor._\n")
    for cand in watch:
        lines.append(
            f"- **{cand['lineage']}** (parent clade {cand['parent_clade']}, "
            f"{len(cand['spike_mutations'])} spike mut) — "
            f"global {pct(cand['frequency']['build'].get('global'))}, "
            f"NA {pct(cand['frequency']['build'].get('north-america'))}, "
            f"EU {pct(cand['frequency']['build'].get('europe'))}; "
            f"MLR ga {cand['mlr_ga'] if cand['mlr_ga'] is not None else 'n/a'}"
        )

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _render_candidate(lines, cand, suggested, header="## "):
    lin = cand["lineage"]
    lines.append(f"\n{header}{lin}  →  suggested name **{suggested}** "
                 f"(parent clade {cand['parent_clade']})")
    if cand.get("alternatives"):
        lines.append(f"\n> Alternative demarcation of the same sweep — choose ONE of "
                     f"{lin} or {', '.join(cand['alternatives'])} to carry {suggested}.")
    if cand["is_recombinant"]:
        lines.append("\n> ⚠ Recombinant (X*) lineage — sanity-check the breakpoint and that the "
                     "defining mutations are inherited cleanly, not split across parents.")
    b = cand["frequency"]["build"]
    lines.append("\n**Frequency (Nextstrain open 6m builds — headline):**")
    lines.append(f"- global {pct(b.get('global'))} | North America {pct(b.get('north-america'))} "
                 f"| Europe {pct(b.get('europe'))}")
    if cand["frequency"]["lapis"]:
        parts = []
        for rk in ("global", "north-america", "europe"):
            if rk in cand["frequency"]["lapis"]:
                v = cand["frequency"]["lapis"][rk]
                warn = "  ⚠low-n" if v["total"] < LOW_DENOMINATOR else ""
                parts.append(f"{rk} {pct(v['freq'])} ({v['n']}/{v['total']}{warn})")
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
            lines.append("\n**MLR fitness/frequency (open):** " + " | ".join(parts))
    lines.append(f"\n**Spike mutations vs parent clade ({len(cand['spike_mutations'])}):** "
                 + (", ".join(cand["spike_mutations"]) or "none"))
    lines.append(f"\n**Global frequency trend (last 6 pivots):** "
                 f"{cand['slope_global']:+.2f}")
    lines.append(f"\n**First seen (open data):** {cand.get('first_sequence') or 'unknown'}"
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
    parser.add_argument("--global-build", default=os.path.join(REPO_ROOT, "auspice", "ncov_open_global_6m.json"))
    parser.add_argument("--global-tip-frequencies", default=os.path.join(REPO_ROOT, "auspice", "ncov_open_global_6m_tip-frequencies.json"))
    parser.add_argument("--clades", default=os.path.join(REPO_ROOT, "defaults", "clades.tsv"))
    parser.add_argument("--display-names", default=os.path.join(REPO_ROOT, "defaults", "clade_display_names.yml"))
    parser.add_argument("--cache-dir", default=os.path.join(REPO_ROOT, "results", "clade_cache"))
    parser.add_argument("--output-md", default=os.path.join(REPO_ROOT, "results", "clade_candidates.md"))
    parser.add_argument("--output-json", default=os.path.join(REPO_ROOT, "results", "clade_candidates.json"))
    parser.add_argument("--date-from", default=None, help="LAPIS window start (default: 180 days before latest pivot)")
    parser.add_argument("--designation-year", type=int, default=None, help="Year for clade names (default: current year, i.e. the designation year)")
    parser.add_argument("--lapis-min", type=float, default=0.15, help="Only run LAPIS for candidates at/above this build frequency")
    parser.add_argument("--skip-lapis", action="store_true", help="Skip exact LAPIS frequency queries")
    parser.add_argument("--skip-regional", action="store_true", help="Skip fetching NA/Europe builds")
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

    print("Loading global open build ...", file=sys.stderr)
    builds = {"global": Build(json.load(open(args.global_build)),
                              json.load(open(args.global_tip_frequencies)))}
    latest_pivot = builds["global"].pivots[-1]

    if not args.skip_regional:
        for region_key in REGIONS:
            try:
                print(f"Fetching {region_key} open build ...", file=sys.stderr)
                tree = fetch_json(BUILD_URL.format(region=region_key), cache_dir, f"{region_key}.json")
                tf = fetch_json(TIP_FREQ_URL.format(region=region_key), cache_dir, f"{region_key}_tf.json")
                builds[region_key] = Build(tree, tf)
            except Exception as exc:
                print(f"  ! could not load {region_key} build ({exc}); "
                      f"regional frequency will fall back to LAPIS", file=sys.stderr)

    mlr = {} if args.skip_mlr else load_mlr(MLR_URL, cache_dir)

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
                              use_lapis=not args.skip_lapis, lapis_min=args.lapis_min,
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

    os.makedirs(os.path.dirname(args.output_md), exist_ok=True)
    write_markdown(args.output_md, flagged, watch, suggested_names,
                   {"ref_updated": ref_updated, "latest_pivot": round(latest_pivot, 3)})
    with open(args.output_json, "w") as fh:
        json.dump({
            "meta": {"ref_updated": ref_updated, "latest_pivot": latest_pivot,
                     "date_from": date_from, "thresholds": {
                         "global": GLOBAL_THRESHOLD, "regional": REGIONAL_THRESHOLD}},
            "flagged": flagged,
            "watch": watch,
        }, fh, indent=2)

    print(f"\nFlagged {len(flagged)} candidate(s); {len(watch)} on the watch list.", file=sys.stderr)
    print(f"Wrote {args.output_md}", file=sys.stderr)
    print(f"Wrote {args.output_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
