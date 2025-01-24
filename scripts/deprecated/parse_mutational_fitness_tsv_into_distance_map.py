"""
Supplementary data S2 from Obermeyer et al (https://www.science.org/doi/10.1126/science.abm1208)
is a TSV table that maps mutations such as "S:D614G" to an estimate of "Δ log R" and is available via GitHub.
Here, we convert this TSV table into a JSON compatable with the augur distance command.
To update model run:

python scripts/developer_scripts/parse_mutational_fitness_tsv_into_distance_map.py

and the version the resulting changes to defaults/mutational_fitness_distance_map.json

Updated model outputs are available at https://github.com/bkotzen/sars-cov2-modeling following:

https://raw.githubusercontent.com/bkotzen/sars-cov2-modeling/main/2024-07-22/PyR0/mutations.tsv

--------------------------------------------------------------

This analysis was removed from the workflow on 2025-01-23
This was drawn from results at https://github.com/bkotzen/sars-cov2-modeling
But this repo hasn't been updated since 2024-07-22
If these results become updated more frequently, we should restore this analysis

This was used in the workflow following:

rule mutational_fitness:
    input:
        tree = "results/{build_name}/tree.nwk",
        alignments = lambda w: rules.translate.output.translations,
        distance_map = config["files"]["mutational_fitness_distance_map"]
    output:
        node_data = "results/{build_name}/mutational_fitness.json"
    benchmark:
        "benchmarks/mutational_fitness_{build_name}.txt"
    log:
        "logs/mutational_fitness_{build_name}.txt"
    params:
        genes = ' '.join(config.get('genes', ['S'])),
        compare_to = "root",
        attribute_name = "mutational_fitness"
    conda:
        config["conda_environment"],
    resources:
        mem_mb=2000
    shell:
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.compare_to} \
            --attribute-name {params.attribute_name} \
            --map {input.distance_map} \
            --output {output} 2>&1 | tee {log}
"""

import argparse
import pandas as pd
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert mutational fitness values to an Augur distance map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", default="https://raw.githubusercontent.com/bkotzen/sars-cov2-modeling/main/2024-07-22/PyR0/mutations.tsv", help="TSV file of mutational effects")
    parser.add_argument("--output", default="defaults/mutational_fitness_distance_map.json", help="JSON file for augur distance")
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
        if "STOP" not in aa_change:
            from_aa = aa_change[0]
            to_aa = aa_change[-1]
            pos_aa = aa_change[1:-1]
            if gene not in structured_mapping:
                structured_mapping[gene] = {}
            if pos_aa not in structured_mapping[gene]:
                structured_mapping[gene][pos_aa] = []
            entry = {"from": from_aa, "to": to_aa, "weight": round(delta_log_R, 10)}
            structured_mapping[gene][pos_aa].append(entry)

    # output this mapping as an augur distance compatable JSON
    # include very slightly negative default to prevent heavily diverged artifactual genomes from
    # appearing as high fitness
    json_output = {"default": -0.003}
    json_output["map"] = structured_mapping

    print("writing mutational_fitness_distance_map.json to defaults/")
    with open(args.output, 'w') as f:
        json.dump(json_output, f, indent=2)
