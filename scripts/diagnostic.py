"""
Use nextclade QC to produce a list of sequences to be excluded.
"""
import argparse
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def datestr_to_ordinal(x, minus_weeks=0):
    try:
        return (datetime.strptime(x,"%Y-%m-%d") - timedelta(weeks=minus_weeks)).toordinal()
    except:
        return np.nan

def earliest_clade_date(Nextstrain_clade, clade_emergence_dates_filename, window_weeks=2):
    clade_dates = pd.read_csv(clade_emergence_dates_filename, index_col="Nextstrain_clade", sep='\t')
    try:
        return datestr_to_ordinal(clade_dates.loc[Nextstrain_clade]['first_sequence'], minus_weeks=window_weeks)
    except:
        return np.nan

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type=str, required=True, help="metadata")
    parser.add_argument("--clade_emergence_dates", type=str, default="defaults/clade_emergence_dates.tsv", help="tsv file with two columns: Nextstrain_clade name and first known sequence for that clade.")
    parser.add_argument("--clock-filter-recent", type=float, default=20, help="maximal allowed deviation from the molecular clock")
    parser.add_argument("--clock-filter", type=float, default=15, help="maximal allowed deviation from the molecular clock")
    parser.add_argument("--snp-clusters", type=int, default=1, help="maximal allowed SNP clusters (as defined by nextclade)")
    parser.add_argument("--rare-mutations", type=int, default=15, help="maximal allowed rare mutations as defined by nextclade")
    parser.add_argument("--clock-plus-rare", type=int, default=25, help="maximal allowed rare mutations as defined by nextclade")
    parser.add_argument("--clade-emergence-window", type=int, default=2, help="number of weeks before official emergence of clade at which sequences can safely be excluded")
    parser.add_argument("--output-exclusion-list", type=str, required=True, help="Output to-be-reviewed addition to exclude.txt")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')

    check_recency = "date_submitted" in metadata
    if check_recency:
        recency_cutoff = (datetime.today() - timedelta(weeks=4)).toordinal()
        recent_sequences = metadata.date_submitted.apply(lambda x: datestr_to_ordinal(x)>recency_cutoff)
    else:
        print("Skipping QC steps which rely on submission recency, as metadata is missing 'date_submitted'")

    check_clade_dates = "Nextstrain_clade" in metadata
    #auto exclude sequences N weeks before their clade emergence
    if check_clade_dates:
        dates = metadata.date.apply(lambda x: datestr_to_ordinal(x))
        clade_dates = metadata.Nextstrain_clade.apply(lambda x: earliest_clade_date(x, args.clade_emergence_dates, window_weeks=args.clade_emergence_window))
    else:
        print("Skipping QC steps which rely on clade-date combinations, as metadata is missing 'Nextstrain_clade'")

    if "clock_deviation" in metadata.columns:
        clock_deviation = np.array([float(x) if isfloat(x) else np.nan for x in metadata.clock_deviation])
    else:
        clock_deviation = np.zeros(len(metadata), dtype=bool)

    if "rare_mutations" in metadata.columns:
        rare_mutations = np.array([float(x) if isfloat(x) else np.nan for x in metadata.rare_mutations])
    else:
        rare_mutations = np.zeros(len(metadata), dtype=bool)

    if "snp_clusters" in metadata.columns:
        snp_clusters = np.array([float(x) if isfloat(x) else np.nan for x in metadata.snp_clusters])
    else:
        snp_clusters = np.zeros(len(metadata), dtype=bool)

    to_exclude = np.zeros_like(clock_deviation, dtype=bool)
    to_exclude |= np.abs(clock_deviation)>args.clock_filter_recent
    if check_recency:
        to_exclude |= (np.abs(clock_deviation)>args.clock_filter)&(~recent_sequences)
        to_exclude |= (rare_mutations>args.rare_mutations)&(~recent_sequences)
        to_exclude |= (np.abs(clock_deviation+rare_mutations)>args.clock_plus_rare)&(~recent_sequences)
    else:
        to_exclude |= (rare_mutations>args.rare_mutations)
        to_exclude |= (np.abs(clock_deviation+rare_mutations)>args.clock_plus_rare)
    if check_clade_dates:
        to_exclude |= dates<clade_dates
    to_exclude |= snp_clusters>args.snp_clusters


    if "QC_mixed_sites" in metadata.columns:
        to_exclude |= metadata.QC_mixed_sites=='bad'

    # write out file with sequences flagged for exclusion
    with open(args.output_exclusion_list, 'w') as excl:
        for s in metadata.loc[to_exclude,'strain']:
            excl.write(f'{s}\n')
