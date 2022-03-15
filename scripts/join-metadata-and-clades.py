#!/usr/bin/env python3
"""
Copied from https://github.com/nextstrain/ncov-ingest/blob/master/bin/join-metadata-and-clades
"""
import argparse
import sys
from datetime import datetime
import pandas as pd
import numpy as np

INSERT_BEFORE_THIS_COLUMN = "pango_lineage"
METADATA_JOIN_COLUMN_NAME = 'strain'
NEXTCLADE_JOIN_COLUMN_NAME = 'seqName'
VALUE_MISSING_DATA = '?'

rate_per_day = 0.0007 * 29903 / 365
reference_day = datetime(2020,1,1).toordinal()

column_map = {
    "clade": "Nextstrain_clade",
    "Nextclade_pango": "Nextclade_pango",
    "totalMissing": "missing_data",
    "totalSubstitutions": "divergence",
    "totalNonACGTNs": "nonACGTN",
    "privateNucMutations.totalUnlabeledSubstitutions":  "rare_mutations",
    "privateNucMutations.totalReversionSubstitutions": "reversion_mutations",
    "privateNucMutations.totalLabeledSubstitutions": "potential_contaminants",
    "qc.overallScore": "QC_overall_score",
    "qc.overallStatus": "QC_overall_status",
    "qc.missingData.status": "QC_missing_data",
    "qc.mixedSites.status": "QC_mixed_sites",
    "qc.privateMutations.status": "QC_rare_mutations",
    "qc.snpClusters.status": "QC_snp_clusters",
    "qc.frameShifts.status": "QC_frame_shifts",
    "qc.stopCodons.status": "QC_stop_codons",
    "frameShifts": "frame_shifts",
    "deletions": "deletions",
    "insertions": "insertions",
    "substitutions": "substitutions",
    "aaSubstitutions": "aaSubstitutions"
}

preferred_types = {
    "divergence": "int32",
    "nonACGTN": "int32",
    "missing_data": "int32",
    "snp_clusters": "int32",
    "rare_mutations": "int32"
}

def reorder_columns(result: pd.DataFrame):
    """
    Moves the new clade column after a specified column
    """
    columns = list(result.columns)
    columns.remove(column_map['clade'])
    insert_at = columns.index(INSERT_BEFORE_THIS_COLUMN)
    columns.insert(insert_at, column_map['clade'])
    return result[columns]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Joins metadata file with Nextclade clade output",
    )
    parser.add_argument("first_file")
    parser.add_argument("second_file")
    parser.add_argument("-o", default=sys.stdout)
    return parser.parse_args()

def datestr_to_ordinal(x):
    try:
        return datetime.strptime(x,"%Y-%m-%d").toordinal()
    except:
        return np.nan

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def main():
    args = parse_args()

    metadata = pd.read_csv(args.first_file, index_col=METADATA_JOIN_COLUMN_NAME,
                           sep='\t', low_memory=False)

    # Check for existing annotations in the given metadata. Skip join with
    # Nextclade QC file, if those annotations already exist and none of the
    # columns have empty values. In the case where metadata were combined from
    # different sources with and without annotations, the "clock_deviation"
    # column will exist but some values will be missing. We handle this case as
    # if the annotations do not exist at all and reannotate all columns. We
    # cannot look for missing values across all expected columns as evidence of
    # incomplete annotations, since a complete annotation by Nextclade will
    # include missing values for some columns by design.
    expected_columns = list(column_map.values()) + ["clock_deviation"]
    existing_annotation_columns = metadata.columns.intersection(expected_columns)
    if len(existing_annotation_columns) == len(expected_columns):
        if metadata["clock_deviation"].isnull().sum() == 0:
            print(f"Metadata file '{args.first_file}' has already been annotated with Nextclade QC columns. Skipping re-annotation.")
            metadata.to_csv(args.o, sep="\t")
            return

    # Read and rename clade column to be more descriptive
    clades = pd.read_csv(args.second_file, index_col=NEXTCLADE_JOIN_COLUMN_NAME,
                         sep='\t', low_memory=False, na_filter = False) \
            .rename(columns=column_map)

    clade_columns = clades.columns.intersection(list(column_map.values()))
    clades = clades[clade_columns]

    # Concatenate on columns
    result = pd.merge(
        metadata, clades,
        left_index=True,
        right_index=True,
        how='left',
        suffixes=["_original", ""],
    )
    all_clades = result.Nextstrain_clade.unique()
    t = result["date"].apply(datestr_to_ordinal)
    div_array = np.array([float(x) if isfloat(x) else np.nan for x in result.divergence])
    offset_by_clade = {}
    for clade in all_clades:
        ind = result.Nextstrain_clade==clade
        if ind.sum()>100:
            deviation = div_array[ind] - (t[ind] - reference_day)*rate_per_day
            offset_by_clade[clade] = np.mean(deviation[~np.isnan(deviation)])

    # extract divergence, time and offset information into vectors or series
    offset = result["Nextstrain_clade"].apply(lambda x: offset_by_clade.get(x, 2.0))
    # calculate divergence
    result["clock_deviation"] = np.array(div_array - ((t-reference_day)*rate_per_day + offset), dtype=int)
    result.loc[np.isnan(div_array)|np.isnan(t), "clock_deviation"] = np.nan

    for col in list(column_map.values()) + ["clock_deviation"]:
        result[col] = result[col].fillna(VALUE_MISSING_DATA)

    # Move the new column so that it's next to other clade columns
    if INSERT_BEFORE_THIS_COLUMN in result.columns:
        result = reorder_columns(result) #.astype(preferred_types)

    result.to_csv(args.o, index_label=METADATA_JOIN_COLUMN_NAME, sep='\t')


if __name__ == '__main__':
    main()
