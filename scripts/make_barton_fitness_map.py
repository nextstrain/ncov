"""
Translate pangolineages from CSV -> JSON for node_data
Note: this should arguably live instead as part of `combine_metadata`,
but this gets particularly complex given the new multiple-inputs logic.
So, for now, following the initial suggestion in the issue.
"""

import json
import pandas as pd

if __name__ == '__main__':

    d = pd.read_csv('https://raw.githubusercontent.com/bartonlab/paper-SARS-CoV-2-inference/main/data/selection-g40-1pct-nonsyn-paper.csv')

    fitness = {}
    for ri, row in d.iterrows():
        pos = row["nucleotide number"]
        if pos not in fitness:
            fitness[pos] = []

        fitness[pos].append({'to':row.nucleotide, "from": row["reference nucleotide"], "weight":row["selection coefficient"]})

    with open('barton_map.json', 'w') as fh:
        json.dump({
            "name": "Barton map",
            "default":0.0,
            "map":{
                "nuc": fitness
            }
        }, fh)
