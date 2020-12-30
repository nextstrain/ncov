
ruleorder: polytomy_tree > refine

rule polytomy_tree:
    message:
        """
        Creating a fake tree (all tips in a polyomy) for auspice to focus on the map
        """
    input:
        alignment = rules.combine_samples.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        tree = "results/{build_name}/tree.nwk",
        node_data = "results/{build_name}/branch_lengths.json"
    threads: 1
    conda: config["conda_environment"]
    run:
        from augur.utils import read_metadata, get_numerical_dates
        import json
        from Bio import SeqIO
        strains = list(SeqIO.index(input.alignment, 'fasta').keys())
        metadata, _ = read_metadata(input.metadata)
        dates = get_numerical_dates(metadata, fmt="%Y-%m-%d", min_max_year=None)
        node_data = {'nodes': {n:{'mutation_length': 1, 'num_date': dates[n]} for n in strains}}
        node_data['nodes']['ROOT'] = {'mutation_length': 0, 'num_date': 2019.9}
        print(f"Single tree with polytomy of {len(strains)} strains")
        with open(output.node_data, 'w') as fh:
            json.dump(node_data, fh, indent=2)
        with open(output.tree, 'w') as fh:
            print(f"({','.join(strains)})ROOT;", file=fh)

ruleorder: dont_incorporate_travel_history > incorporate_travel_history

rule dont_incorporate_travel_history:
    input:
        auspice_json = rules.export.output.auspice_json,
    output:
        auspice_json = "results/{build_name}/ncov_with_accessions_and_travel_branches.json"
    shell:
        """
        cp {input.auspice_json} {output.auspice_json}
        """
