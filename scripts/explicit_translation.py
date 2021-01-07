import argparse
import json
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from treetime import TreeAnc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add translations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="input tree")
    parser.add_argument('--translation', type=str, required=True, help="amino acid alignment")
    parser.add_argument('--gene', type=str, required=True, help="amino acid alignment")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')
    leafs = {n.name for n in T.get_terminals()}
    seqs = []
    for s in SeqIO.parse(args.translation, 'fasta'):
        if s.id in leafs:
            seqs.append(s)


    tt = TreeAnc(tree=T, aln=MultipleSeqAlignment(seqs), alphabet='aa')

    tt.infer_ancestral_sequences(reconstruct_tip_states=True)

    node_data = {}

    for n in tt.tree.find_clades():
        node_data[n.name] = {"aa_muts":{args.gene:[f"{a}{p+1}{d}" for a,p,d in n.mutations]}}

    with open(args.output, 'w') as fh:
        json.dump({"nodes":node_data}, fh)
