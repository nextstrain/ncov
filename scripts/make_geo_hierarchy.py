"""
for each country, region, division, location glob the containing region, country
etc out of the metadata. 
"""

import argparse
import json
import pandas as pd
from collections import defaultdict

def get_max(counts):
    items = sorted(counts.items(), key=lambda x:x[1])
    return items[-1][0]


def majority_rule(d):
    max_only = {}
    for k,v in d.items():
        max_only[k] = {h:get_max(counts) for h, counts in v.items()}

    return max_only

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make region hierarchy",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--output", type=str, required=True, help="adjusted metadata")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, delimiter='\t')

    regions = {}
    countries = {}
    divisions = {}
    locations = {}

    for ri, row in metadata.iterrows():
        region, country, division, location = map(str, [row['region'], row['country'], row['division'], row['location']])
        if region!='nan':
            if region not in regions:
                regions[region] = {'region':defaultdict(int)}
            regions[region]['region'][region] += 1

        if region!='nan' and country!='nan':
            if country not in countries:
                countries[country] = {'region':defaultdict(int), 'country':defaultdict(int)}
            countries[country]['region'][region] += 1
            countries[country]['country'][country] += 1

        if region!='nan' and country!='nan' and division!='nan':
            if division not in divisions:
                divisions[division] = {'region':defaultdict(int), 'country':defaultdict(int), 'division':defaultdict(int)}
            divisions[division]['region'][region] += 1
            divisions[division]['country'][country] += 1
            divisions[division]['division'][division] += 1

        if region!='nan' and country!='nan' and division!='nan' and location!='nan':
            if location not in locations:
                locations[location] = {'region':defaultdict(int), 'country':defaultdict(int), 'division':defaultdict(int), 'location':defaultdict(int)}
            locations[location]['region'][region] += 1
            locations[location]['country'][country] += 1
            locations[location]['division'][division] += 1
            locations[location]['location'][location] += 1

    

    geo_hierarchy = {"region": majority_rule(regions), 
                     "country": majority_rule(countries), 
                     "division": majority_rule( divisions), 
                     "location": majority_rule(locations)
                    }

    with open(args.output, 'w') as fh:
        json.dump(geo_hierarchy, fh)
