---
name: designate-clades
description: >-
  Propose and draft new SARS-CoV-2 Nextstrain clade designations from open data.
  TRIGGER when the user wants to find/designate/elevate a new clade, check whether
  any Pango lineage has crossed the clade-naming threshold, refresh clade candidates,
  or draft a clade-designation PR in the nextstrain/ncov repo. Runs
  scripts/propose_clades.py, interprets the candidate dossier, writes the five
  defaults/* edits, validates them, and drafts the PR.
---

# Designating SARS-CoV-2 Nextstrain clades

Clade designation is a judgement call on top of a quantitative threshold. The
`scripts/propose_clades.py` tool does the *evidence gathering*; you make the
*decisions* (which lineage to elevate, where to draw the line, the name, the prose)
and write the edits. The tool never edits `defaults/` — that is your job here.

## Designation criteria

A candidate Pango lineage qualifies when it **both**:
1. carries **>=1 spike (S:) mutation** relative to its parent Nextstrain clade, and
2. has risen to **>30% regional** (North America or Europe) **or >20% global** frequency.

Only ever elevate an **existing Pango lineage** to a clade — keep a strict 1:1
Pango↔Nextstrain-clade mapping. The hard part is usually *which* lineage in a nested
chain (e.g. XFG.1 vs XFG.1.1) to elevate: **prefer the branch with the larger fitness
differential**, visible as a jump in MLR growth advantage and corroborated by the
frequency trajectory.

Open data is sparse and ~2 months lagged in 2026 (GISAID API access was lost — see the
2026-02-24 blog). Treat every percentage as approximate and triangulate the three
frequency reads the dossier gives you (build tip-freq, LAPIS exact, MLR modeled).

## Step 1 — Generate the dossier

```
python3 scripts/propose_clades.py --data-source open               # global + NA + Europe
python3 scripts/propose_clades.py --data-source gisaid --group blab # global + all 6 regions
```

Everything is fetched live (cached under `results/clade_cache/`; `--refresh` to re-pull): the
Nextclade reference tree, the Nextstrain 6m builds, forecasts-ncov MLR fitness, and — for **open**
only — LAPIS exact frequencies. **gisaid** auto-discovers the group's latest dated Groups builds, so
the date is handled for you (`--gisaid-date` pins an older one). open is GISAID-free but
North-America-dominated; gisaid covers all regions and surfaces sweeps open misses.

Frequency is judged conservatively: a fixed-N=50 window of date-sorted tips slides across a
~6mo→3wk band and a candidate flags only if the **one-sided 95% Wilson lower bound** clears the
threshold (>20% global / >30% regional) in some window — so each flag is a *confident* majority
(effective bar ≈46% regional / ≈34% global at N=50), reported as `peak X% [n, as of <date>, lo Y%] ·
now Z%`. A geography with <50 band tips reads `insufficient`. So a flag already means the data backs
it; your job is the demarcation/naming judgement, not re-litigating noise.

Outputs:
- `results/clade_candidates.md` — human-readable dossier (read this)
- `results/clade_candidates.json` — same data, structured (paste-ready edits live here)

## Step 2 — Interpret each flagged candidate

For every candidate under "Flagged", read the evidence block and cross-check against the
live resources before trusting the numbers:
- **Nextstrain builds:** https://nextstrain.org/ncov/open/global/6m ,
  `.../north-america/6m` , `.../europe/6m`
- **Cov-Spectrum** (exact, per lineage): e.g.
  `https://cov-spectrum.org/explore/World/AllSamples/Past6M/variants?nextcladePangoLineage=<LIN>*`
- **MLR fitness/forecast:** https://nextstrain.github.io/forecasts-ncov/
- **Nextclade reference tree** (mutations): https://nextstrain.org/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs

Use WebFetch on the LAPIS open API for an exact recent read if you want a tighter window
than the dossier's 180 days, e.g.
`https://lapis.cov-spectrum.org/open/v2/sample/aggregated?nextcladePangoLineage=<LIN>*&region=North%20America&dateFrom=<YYYY-MM-DD>`
(divide by the same query without the lineage filter for the denominator).

Weigh, since a flag already passed the windowed CI:
- **only one small region clears** (e.g. a 56% Oceania-only majority) — decide whether it warrants a
  *global* clade name or is just locally dominant; the per-geography lines and `now` values show reach;
- the **peak is months old and `now` has collapsed** — still designation-worthy retrospectively, but
  note it in the PR;
- the lineage is **recombinant (X*)** — confirm the breakpoint and that the defining mutations are
  inherited cleanly;
- for **open**, LAPIS is shown as a secondary exact read; large build-vs-LAPIS gaps are expected
  (LAPIS is raw/country-biased, the build is population-weighted) — prefer the build for the call.

## Step 3 — Decide the demarcation

Look at the "Where to draw the line" table. Walk the nested chain (parent → candidate →
children) and pick the lineage where the **MLR growth advantage jumps most** (largest `Δga`
/ highest `ga`) and whose frequency is actually rising (`trend` > 0). That is the branch
carrying the real fitness differential. Elevating a parent that is mostly one fast child, or
a child that is only a slice of an already-rising parent, are both mistakes — the table is
there to tell those apart.

If the answer is "this is just slicing the current dominant clade finer with no real fitness
story," it is fine to **not** designate. Note it on the watch list and move on.

## Step 4 — Write the five `defaults/` edits

The dossier prints a paste-ready block per candidate. Apply it, then **trim the
`clades.tsv` nuc rows to the characteristic subset** (existing clades use ~2–4 defining
mutations, not the full accumulated set — drop homoplasic / unstable sites). Confirm the
chosen mutations are present in essentially all members of the lineage.

1. `defaults/clades.tsv` — `<NAME>\tclade\t<PARENT>` + a few `<NAME>\tnuc\t<site>\t<alt>` rows
2. `defaults/clade_display_names.yml` — `<NAME>: <NAME> (<Pango>)`
3. `defaults/clade_hierarchy.tsv` — `<NAME>\t<PARENT>\t<WHO>`
4. `defaults/clade_emergence_dates.tsv` — `<NAME>\t<YYYY-MM-01>` (first-seen from the dossier)
5. `defaults/color_ordering.tsv` — `clade_membership\t<NAME> (<Pango>)`

The suggested `<NAME>` (e.g. `26A`) uses the **designation year** — the year the clade is being
named, i.e. the current year — *not* the lineage's emergence year. It is the first unused letter
for that year (so the first clade named in 2026 is `26A`, the next `26B`, and so on). The
emergence date still goes in `clade_emergence_dates.tsv` as the lineage's first-seen date. Only
touch `emerging_lineages.tsv` for recombinant/special lineages, removing a lineage there if you
are promoting it to a full clade (see PR #1152).

## Step 5 — Validate

Confirm the new definition parses and inherits correctly through augur:

```
augur shell .  # or use the repo's nextstrain runtime
scripts/expand-clade-definitions defaults/clades.tsv | grep "^<NAME>"
```

(`scripts/expand-clade-definitions` imports `augur.clades.read_in_clade_definitions`; run it
with a Python that has augur installed.) Ideally also re-run the `clades` rule against the
local global build and confirm the expected tips pick up `<NAME>` and it renders in Auspice.

## Step 6 — Draft the PR

Open a PR in the house style of #1149 / #1152 / #1158 / #1192. The justification should state:
- the new clade name and its Pango lineage, and the parent clade it descends from;
- the **spike mutations** that define it;
- the **frequency** evidence with regions and numbers (which criterion arm it meets);
- the **MLR fitness** (growth advantage) and that it is rising;
- a one-line note on the demarcation choice (why this lineage and not its parent/child).

Keep the dossier's exact numbers in the PR body so reviewers can trace them. Do not push or
open the PR without the user's go-ahead.
