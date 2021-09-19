"""
Supplementary data S2 from Obermeyer et al (https://www.medrxiv.org/content/10.1101/2021.09.07.21263228v1)
is a TSV table that maps mutations such as "S:D614G" to an estimate of "Δ log R".
Here, we convert this TSV table into a JSON compatable with the augur distance command.
"""
import argparse
import pandas as pd
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="TSV file of mutational effects")
    parser.add_argument("--output", required=True, help="JSON file for augur distance")
    args = parser.parse_args()

    # collect simple string mapping from TSV, ie
    # 'S:A522V': -0.00661378
    # 'ORF1a:T4304I': 0.00353199
    string_mapping = {}
    if args.input:
        df = pd.read_csv(args.input, delimiter='\t')
        for index, row in df.iterrows():
            string_mapping[row["mutation"]] = float(row["Δ log R"])

    # convert simple string mapping into structured mapping required by augur distance, ie
    # "map": {
    #   "S": {
    #     "522": [
    #       {
    #         "from": "A",
    #         "to": "V",
    #         "weight": -0.00661378
    #       },
    structured_mapping = {}
    for mutation, delta_log_R in string_mapping.items():
        gene, aa_change = mutation.split(":")
        from_aa = aa_change[0]
        to_aa = aa_change[-1]
        pos_aa = aa_change[1:-1]
        if gene not in structured_mapping:
            structured_mapping[gene] = {}
        if pos_aa not in structured_mapping[gene]:
            structured_mapping[gene][pos_aa] = []
        entry = { "from": from_aa, "to": to_aa, "weight": delta_log_R }
        structured_mapping[gene][pos_aa].append(entry)

    # output this mapping as an augur distance compatable JSON
    # include very slightly negative default to prevent heavily diverged artifactual genomes from
    # appearing as high fitness
    json_output = { "default": -0.003 }
    json_output["map"] = structured_mapping

    with open(args.output, 'w') as f:
        json.dump(json_output, f, indent=2)
