import argparse
import json
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from treetime import TreeAnc
from augur.utils import load_features


def annotation_json(features, reference):
    annotations = {}
    for fname, feat in features.items():
        annotations[fname] = {'seqid':reference.id,
                              'type':feat.type,
                              'start':int(feat.location.start)+1,
                              'end':int(feat.location.end),
                              'strand': '+' if feat.location.strand else '-'}
    annotations['nuc'] = {'seqid':reference.id,
                            'type':'source',
                            'start': 1,
                            'end': len(reference),
                            'strand': '+'}
    return annotations


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add translations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="input tree")
    parser.add_argument('--reference', type=str, required=True, help="reference genbank sequence")
    parser.add_argument('--translations', type=str,  nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--genes', type=str, nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    genes = args.genes if type(args.genes)==list else [args.genes]
    translations = args.translations if type(args.translations)==list else [args.translations]
    ref = SeqIO.read(args.reference, format='genbank')
    features = load_features(args.reference)

    if not set(features.keys())==set(args.genes):
        print("WARNING: supplied genes don't match the annotation")
        print("the following features are in the annotation by not supplied as genes:", set(features.keys()).difference(args.genes))
        print("the following features are in the supplied as genes but not the annotation:", set(args.genes).difference(features.keys()))

    T = Phylo.read(args.tree, 'newick')
    leafs = {n.name for n in T.get_terminals()}

    node_data = {}
    root_sequence_translations = {}
    for gene, translation in zip(genes, translations):
        seqs = []
        for s in SeqIO.parse(translation, 'fasta'):
            if s.id in leafs:
                seqs.append(s)


        tt = TreeAnc(tree=T, aln=MultipleSeqAlignment(seqs), alphabet='aa')

        tt.infer_ancestral_sequences(reconstruct_tip_states=True)
        root_sequence_translations[gene] = tt.sequence(tt.tree.root, as_string=True, reconstructed=True)

        with open(translation.replace('.fasta', '_withInternalNodes.fasta'), 'w') as fh:
            for n in tt.tree.find_clades():
                if n.name not in node_data:
                    node_data[n.name] = {"aa_muts":{}}
                if len(n.mutations):
                    node_data[n.name]["aa_muts"][gene] = [f"{a}{p+1}{d}" for a,p,d in n.mutations]
                fh.write(f">{n.name}\n{tt.sequence(n, as_string=True, reconstructed=True)}\n")

    annotations = annotation_json(features, ref)
    with open(args.output, 'w') as fh:
        json.dump({"nodes":node_data, "annotations":annotations, "reference":root_sequence_translations}, fh)
