"""
Add column to metadata to denote 'focal' samples based on supplied region
Rewrite location, division and country for non-focal samples to be region
Rewrite division_exposure and country_exposure for non-focal samples to be region_exposure
"""

import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add column to metadata to denote 'focal' samples based on supplied region",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--country", type=str, required=True, help="focal country")
    parser.add_argument("--output", type=str, required=True, help="adjusted metadata")
    args = parser.parse_args()

    fix_casing = {
        "iceland": "Iceland",
        "new zealand": "New Zealand"
        # "africa": "Africa",
        # "europe": "Europe",
        # "north america": "North America",
        # "oceania": "Oceania",
        # "south america": "South America"
    }
    focal_country = fix_casing[args.country]


    print("Adjusting metadata for focal country", args.country)

    metadata = pd.read_csv(args.metadata, delimiter='\t')

    metadata.loc[metadata.country != focal_country, 'location'] = ""
    metadata.loc[metadata.country != focal_country, 'division'] = metadata.country
    metadata.loc[metadata.country != focal_country, 'division_exposure'] = metadata.country_exposure
    metadata.loc[(metadata.country == focal_country) & (metadata.country_exposure != focal_country), 'division_exposure'] = metadata.country_exposure

    metadata.to_csv(args.output, index=False, sep="\t")
