import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Check for missing colors & locations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, nargs='+', required=True, help="input region adjusted metadata")
    parser.add_argument('--colors', type=str, nargs='+', required=True, help="input region specific color file")
    parser.add_argument('--latlong', type=str, required=True, help="input lat-long file")
    args = parser.parse_args()

    things_to_exclude_orig = ['Africa', 'Asia', 'South America', 'Europe', 
                              'North America', 'Oceania', 'Grand princess cruise ship',
                              'diamond princess']
    things_to_exclude = [x.lower() for x in things_to_exclude_orig]

    all_metadatas = [pd.read_csv(met, delimiter='\t') for met in args.metadata]
    metadata = pd.concat(all_metadatas, sort=False)
    all_colors = [pd.read_csv(col, delimiter='\t', header=None) for col in args.colors]
    colors = pd.concat(all_colors, sort=False)

    latlong = pd.read_csv(args.latlong, delimiter='\t', header=None)

    for geo_value in ['location', 'division', 'country']:
        locs_w_color_orig = colors.loc[colors[0]==geo_value,1].values
        locs_w_color = [x.lower() for x in locs_w_color_orig]
        locs_w_latlong_orig = latlong.loc[latlong[0]==geo_value,1].values
        locs_w_latlong = [x.lower() for x in locs_w_latlong_orig]
        locs_in_meta_orig = [x for x in metadata[geo_value].unique() if not pd.isna(x)]
        locs_in_meta = [x.lower() for x in locs_in_meta_orig]

        missing_color_locs = [loc for loc in locs_in_meta if loc not in locs_w_color and loc not in things_to_exclude]
        if missing_color_locs:
            print("The following {} are missing colors:".format(geo_value))
            print(missing_color_locs)
            print("\n")

        if geo_value != 'country':
            missing_latlong_locs = [loc for loc in locs_in_meta if loc not in locs_w_latlong and loc not in things_to_exclude]
            if missing_latlong_locs:
                print("The following {} are missing lat-long values:".format(geo_value))
                print(missing_latlong_locs)
                print("\n")
    
    print("Please remember this does *not* check lat-longs for countries!!")
