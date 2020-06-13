"""
Add column to metadata to denote 'focal' samples based on supplied geo resolution
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
    parser.add_argument("--focal-resolution", type=str, required=False, help="focal geo resolution, eg region, country, division")
    parser.add_argument("--focal-label", type=str, required=False, help="focal geo label, eg North America, USA, Washington")
    parser.add_argument("--output", type=str, required=True, help="adjusted metadata")
    args = parser.parse_args()

    focal_label = args.focal_label
    focal_resolution = args.focal_resolution
    focal_resolution_exposure = args.focal_resolution + "_exposure"

    print("Adjusting metadata for focal geography", focal_label, "with resolution", focal_resolution)

    metadata = pd.read_csv(args.metadata, delimiter='\t')
    metadata.insert(12, 'focal', True)

    metadata.loc[metadata[focal_resolution] != focal_label, 'focal'] = False
    metadata.loc[metadata[focal_resolution] != focal_label, 'location'] = ""
    metadata.loc[metadata[focal_resolution] != focal_label, 'division'] = metadata.region
    metadata.loc[metadata[focal_resolution] != focal_label, 'country'] = metadata.region
    metadata.loc[metadata[focal_resolution] != focal_label, 'division_exposure'] = metadata.region_exposure
    metadata.loc[metadata[focal_resolution] != focal_label, 'country_exposure'] = metadata.region_exposure
    metadata.loc[(metadata[focal_resolution] == focal_label) & (metadata[focal_resolution_exposure] != focal_label), 'division_exposure'] = metadata.region_exposure
    metadata.loc[(metadata[focal_resolution] == focal_label) & (metadata[focal_resolution_exposure] != focal_label), 'country_exposure'] = metadata.region_exposure
    metadata.loc[(metadata[focal_resolution] == focal_label) & (metadata.division_exposure.isna()), 'division_exposure'] = metadata.division
    metadata.loc[(metadata[focal_resolution] == focal_label) & (metadata.country_exposure.isna()), 'country_exposure'] = metadata.country

    metadata.to_csv(args.output, index=False, sep="\t")
