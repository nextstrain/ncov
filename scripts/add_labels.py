import argparse
import json
from Bio import Phylo
from collections import defaultdict

def attach_labels(d, labeled_nodes):
    if "children" in d:
        for c in d["children"]:
            if c["name"] in labeled_nodes:
                if "labels" not in c["branch_attrs"]:
                    c["branch_attrs"]["labels"] = {}
                c['branch_attrs']['labels']['mlabel'] = labeled_nodes[c["name"]][0]
                print(c['branch_attrs']['labels'])
            attach_labels(c, labeled_nodes)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Remove extraneous colorings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--tree', type=str, required=True, help="tree file")
    parser.add_argument('--clades', type=str, required=True, help="clades")
    parser.add_argument('--mutations', type=str, required=True, help="mutations")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')

    with open(args.mutations, "r") as f:
        mutation_json = json.load(f)['nodes']

    with open(args.clades, "r") as f:
        clades_json = json.load(f)['nodes']

    with open(args.input, "r") as f:
        input_json = json.load(f)

    nodes = {}
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            n.tip_count=1
        else:
            n.tip_count = sum([c.tip_count for c in n])
        nodes[n.name] = {'tip_count':n.tip_count}

    labels = defaultdict(list)
    for node in nodes:
        for m in mutation_json[node]['muts']:
            if m[0] in 'ACGT' and m[-1] in 'ACGT':
                clade = clades_json[node]['clade_membership']
                tmp_label = (clade, m)
                labels[tmp_label].append((node, nodes[node]['tip_count']))

    labeled_nodes = defaultdict(list)
    for label in labels:
        node = sorted(labels[label], key=lambda x:-x[1])[0]
        labeled_nodes[node[0]].append('/'.join(label))

    attach_labels(input_json["tree"], labeled_nodes)

    with open(args.output, 'w') as f:
        json.dump(input_json, f, indent=2)
