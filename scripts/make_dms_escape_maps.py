import pandas as pd
from collections import defaultdict
import json
import numpy as np


data = pd.read_csv(r"https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_data.csv")

maps = {x:defaultdict(list) for x in data.condition_subtype.unique()}
sum_maps = {x:defaultdict(lambda: defaultdict(float)) for x in data.condition_subtype.unique()}
mean_maps = {x:defaultdict(lambda: defaultdict(float)) for x in data.condition_subtype.unique()}

for c, row in data.iterrows():
    maps[row.condition_subtype][(row.site, row.wildtype, row.mutation)].append(row.mut_escape)
    sum_maps[row.condition_subtype][row.site][row.condition] = row.normalized_site_total_escape
    mean_maps[row.condition_subtype][row.site][row.condition] = row.normalized_site_mean_escape

average_maps = {x:{} for x in maps}
for x in maps:
    for pos, wt, mut in maps[x]:
        val = np.mean(maps[x][(pos, wt, mut)])
        if val>0.05:
            if pos not in average_maps[x]:
                average_maps[x][pos]={}
            average_maps[x][pos][(wt, mut)]=val

for x in average_maps:
    out_json = {"default":0,
                "map":{"S":{ pos:[{'from':k[0], 'to':k[1], 'weight':y} for k,y in average_maps[x][pos].items()]
                            for pos in average_maps[x]}},
                "name":x}

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)

for x in sum_maps:
    out_json = {"default":0,
                "map":{"S":{}},
                "name":x}
    for pos in sum_maps[x]:
        out_json["map"]["S"][pos] = np.mean(list(sum_maps[x][pos].values()))

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_total_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)


for x in mean_maps:
    out_json = {"default":0,
                "map":{"S":{}},
                "name":x}
    for pos in mean_maps[x]:
        out_json["map"]["S"][pos] = np.mean(list(mean_maps[x][pos].values()))

    with open(f"defaults/distance_maps/{x.replace(' ', '_')}_mean_dms.json", 'w') as fh:
        json.dump(out_json, fh, indent=2)

