from Bio import SeqIO, Phylo
import json
import yaml
import sys
import argparse

def find_matching_nodes(clades_fname, basal_clade_label, tree_fname):
    basal_node_name = None
    with open(clades_fname) as fh:
        for name, node_data in json.load(fh)['nodes'].items():
            if node_data.get('clade_annotation', '') == basal_clade_label:
                basal_node_name = name
                break
    if not basal_node_name:
        print(f"WARNING: no node found with a clade_annotation of {basal_clade_label}. This script will proceed, but no levels will be exported.")
        return set()
    print(f"Node representing {basal_clade_label}: {basal_node_name}")
    T = Phylo.read(tree_fname, 'newick')
    basal_node = T.find_any({"name": basal_node_name})
    if not basal_node:
        print(f"ERROR: {basal_node_name} not found in provided tree") # this should be fatal as it indicates a mismatch of provided inputs
        sys.exit(2)
    names = set([basal_node_name]) # include parent (the basal_clade_label defining branch)
    for n in basal_node.find_clades():
        names.add(n.name)

    return names

def classify_into_levels(spike_seq, rbd_mutations):
    level_num = 0
    codons = []
    calls = []
    for idx in range(0, len(rbd_mutations)):
        (basal_codon, position, alts) = rbd_mutations[idx]
        aa = spike_seq[position-1] # humans use 1-based positions, python uses 0-based
        codons.append(aa)
        if aa==basal_codon:
            calls.append('basal')
        elif aa in alts:
            calls.append('alts')
            level_num +=1
        elif aa=='X' or aa=='-':
            calls.append('X-')
        else:
            calls.append('other')
    return (level_num, codons, calls)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign (omicron) levels to strains",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--config', type=str, metavar="YAML", required=True, help="config defining the RBD mutations")
    parser.add_argument('--spike-alignment', type=str, metavar="FASTA", required=True, help="input spike alignment (usually from nextclade)")
    parser.add_argument('--output-node-data', type=str, metavar="JSON", required=True, help="output node-data JSON")
    parser.add_argument('--clades-node-data', type=str, metavar="JSON", required=False, help="Clade definitions JSON. Only descendants of the basal-clade-label will be labelled.")
    parser.add_argument('--basal-clade-label', type=str, metavar="STRING", required=False, help="Used in conjunction with --clades-node-data and --tree")
    parser.add_argument('--tree', type=str, metavar="NEWICK", required=False, help="Used in conjunction with --clades-node-data and --basal-clade-label")
    args = parser.parse_args()

    if args.clades_node_data and args.basal_clade_label and args.tree:
        selected_nodes = find_matching_nodes(args.clades_node_data, args.basal_clade_label, args.tree)
        use_node = lambda name: name in selected_nodes
    elif args.clades_node_data or args.tree or args.basal_clade_label:
        print("ERROR: if you provide --clades-node-data you must provide --tree, and vice-versa")
        sys.exit(2)
    else:
        use_node = lambda name: True

    with open(args.config, "r") as stream:
        config = yaml.safe_load(stream)
    try:
        rbd_mutations=config['rbd_mutations']
    except KeyError:
        print("ERROR: the provided yaml config did not define a `rbd_mutations` list")
        sys.exit(2)

    spike_aln = {}
    with open(args.spike_alignment) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if use_node(record.id):
                spike_aln[record.id] = record.seq

    node_data = {
        "nodes": {}, # encode the levels for augur export
        "rbd_level_details": {} # for more info, as needed
    }

    for name, seq in spike_aln.items():
        (level_num, codons, calls) = classify_into_levels(seq, rbd_mutations)
        node_data['nodes'][name] = {'rbd_level': level_num}
        node_data['rbd_level_details'][name] = ", ".join([f"S:{x[0][1]}{x[1]} ({x[2]})" for x in zip(rbd_mutations, codons, calls)])

    with open(args.output_node_data, 'w') as fh:
        json.dump(node_data, fh, indent=2)