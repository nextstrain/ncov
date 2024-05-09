import argparse
import itertools
import numpy as np
import pandas as pd
from augur.dates import get_iso_year_week


def export_weights(output):
    # Read data from preprocessed OWID case counts
    # <https://github.com/nextstrain/forecasts-ncov/blob/-/ingest/bin/fetch-ncov-global-case-counts>
    df = pd.read_csv("https://data.nextstrain.org/files/workflows/forecasts-ncov/cases/global.tsv.gz", sep='\t')

    # Rename columns to match names in metadata
    column_name_map = {
        'location': 'country',
        'cases': 'weight',
    }
    df = df.rename(columns=column_name_map)

    # Add groups that are missing since case counts of 0 are filtered out in the
    # forecasts-ncov script. Some may be truly zero while others may be lack of
    # data. These will be filled in further down.
    group_by = ['country', 'date']
    all_dates = df['date'].unique()
    all_countries = df['country'].unique()
    all_groups = set(itertools.product(all_countries, all_dates))
    present_groups = set(df.set_index(group_by).index)
    missing_groups = all_groups - present_groups
    missing_groups_df = pd.DataFrame(missing_groups, columns=group_by)
    missing_groups_df['weight'] = np.nan  # oddly, pd.NA doesn't work here
    df = pd.merge(df, missing_groups_df, how='outer', on=[*group_by, 'weight'])

    # Fill in weights for missing groups. Notes:
    # 1. For missing data flanked by weeks with an available case count, this
    #    uses linear interpolation to "guess" the case count.
    # 2. Countries with abrupt reporting start or cutoff will extrapolate the
    #    first/last available week for the missing weeks. It may be possible to
    #    interpolate from 0 on the very end instead of extrapolating the
    #    constant, but that seems a bit more difficult and not sure if it's any
    #    better.
    # 3. It looks like some case counts don't represent cases within that week
    #    but rather cumulative since the last reported week. Unless that can be
    #    assumed for every week that follows a gap in data, I don't think
    #    there's anything that can be done here.
    df = df.sort_values(group_by)
    df.set_index(group_by, inplace=True)
    df = df.groupby('country').apply(lambda group: group.interpolate(method='linear', limit_direction='both'))
    df.reset_index(inplace=True)

    # Convert YYYY-MM-DD to YYYY-WW
    # Inspired by code in augur.filter.subsample.get_groups_for_subsampling
    # <https://github.com/nextstrain/augur/blob/60a0f3ed2207c5746aa6fc1aa29ab3f75990cb9f/augur/filter/subsample.py#L17>
    temp_prefix = '__ncov_date_'
    temp_date_cols = [f'{temp_prefix}year', f'{temp_prefix}month', f'{temp_prefix}day']
    df_dates = df['date'].str.split('-', n=2, expand=True)
    df_dates = df_dates.set_axis(temp_date_cols[:len(df_dates.columns)], axis=1)

    for col in temp_date_cols:
        df_dates[col] = pd.to_numeric(df_dates[col], errors='coerce').astype(pd.Int64Dtype())

    # Extend metadata with generated date columns
    # Drop the date column since it should not be used for grouping.
    df = pd.concat([df.drop('date', axis=1), df_dates], axis=1)

    df['week'] = df.apply(lambda row: get_iso_year_week(
        row[f'{temp_prefix}year'],
        row[f'{temp_prefix}month'],
        row[f'{temp_prefix}day']
        ), axis=1
    )

    # Output an ordered subset of columns
    df = df[['country', 'week', 'weight']]
    df.to_csv(output, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create weights file",
    )

    parser.add_argument('--output', type=str, metavar="FILE", required=True, help="Path to output weights file")
    args = parser.parse_args()

    export_weights(args.output)
