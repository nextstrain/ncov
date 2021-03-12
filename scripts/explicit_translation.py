
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
    parser.add_argument('--translations', type=str,  nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--genes', type=str, nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    genes = args.genes if type(args.genes)==list else [args.genes]
    translations = args.translations if type(args.translations)==list else [args.translations]

    T = Phylo.read(args.tree, 'newick')
    leafs = {n.name for n in T.get_terminals()}

    node_data = {}
    for gene, translation in zip(genes, translations):
        seqs = []
        for s in SeqIO.parse(translation, 'fasta'):
            if s.id in leafs:
                seqs.append(s)


        tt = TreeAnc(tree=T, aln=MultipleSeqAlignment(seqs), alphabet='aa')

        tt.infer_ancestral_sequences(reconstruct_tip_states=True)

        with open(translation.replace('.fasta', '_withInternalNodes.fasta'), 'w') as fh:
            for n in tt.tree.find_clades():
                if n.name not in node_data:
                    node_data[n.name] = {"aa_muts":{}}
                node_data[n.name]["aa_muts"][gene] = [f"{a}{p+1}{d}" for a,p,d in n.mutations]
                fh.write(f">{n.name}\n{tt.sequence(n, as_string=True, reconstructed=True)}\n")


    with open(args.output, 'w') as fh:
        json.dump({"nodes":node_data}, fh)
