
from __future__ import annotations
import argparse
import json
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from treetime import TreeAnc


def read_gff(fname):
    try:
        from BCBio import GFF #Package name is confusing - tell user exactly what they need!
    except ImportError:
        print("ERROR: Package BCBio.GFF not found! Please install using \'pip install bcbio-gff\' before re-running.")
        return None

    features = {}
    with open(fname, encoding='utf-8') as in_handle:
        for rec in GFF.parse(in_handle):
            for feat in rec.features:
                if "gene_name" in feat.qualifiers:
                    fname = feat.qualifiers["gene_name"][0]
                if fname:
                    features[fname] = feat

    return features

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add translations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True, help="input tree")
    parser.add_argument('--annotation', type=str, required=True, help="gff annotation file")
    parser.add_argument('--reference', type=str, required=True, help="reference fasta sequence")
    parser.add_argument('--translations', type=str,  nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--genes', type=str, nargs='+', required=True, help="amino acid alignment")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    genes = args.genes if type(args.genes)==list else [args.genes]
    translations = args.translations if type(args.translations)==list else [args.translations]
    features = read_gff(args.annotations)
    ref = SeqIO.read(args.reference, format='fasta')

    if not set(features.keys())==set(args.genes):
        print("ERROR: supplied genes don't match the annotation")
        exit(1)

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

    annotation = annotation_json(features, ref)
    with open(args.output, 'w') as fh:
        json.dump({"nodes":node_data, "annotation":annotation}, fh)
