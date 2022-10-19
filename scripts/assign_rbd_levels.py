from Bio import SeqIO
import json
import yaml
import sys
import argparse

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
    args = parser.parse_args()

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