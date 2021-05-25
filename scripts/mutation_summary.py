import argparse, os, glob
from augur.io import open_file
from Bio import SeqIO, SeqFeature, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd



def read_reference(fname, genemap):
    try:
        ref = str(SeqIO.read(fname, 'fasta').seq)
    except:
        with open(fname, 'r') as fh:
            ref = "".join([x.strip() for x in fh])

    translations = {}
    with open(genemap, 'r') as fh:
        for line in fh:
            if line[0]=='#':
                continue
            entries = [x.strip() for x in line.strip().split('\t')]
            start = int(entries[3])
            end = int(entries[4])
            strand = entries[6]
            attributes = {x.split()[0]:' '.join(x.split()[1:]) for x in entries[8].split(';')}
            if 'gene_name' in attributes:
                name = attributes['gene_name'].strip('"')
            else:
                name = None
            translation = Seq.translate(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start-1, end, strand=-1 if strand=='-' else 1)).extract(ref))
            translations[name] = str(translation)

    return {"nuc":ref, "translations":translations}

def get_differences(s1,s2, ambiguous='N'):
    s2_rstrip = s2.rstrip('-')
    s2_lrstrip = s2_rstrip.lstrip('-')
    s2_trimmed = ambiguous*(len(s2_rstrip) - len(s2_lrstrip)) + s2_lrstrip + ambiguous*(len(s2)-len(s2_rstrip))
    return [(a,p+1,d) for p, (a,d) in enumerate(zip(s1,s2_trimmed)) if d not in [a, ambiguous]]

def to_mutations(aln_file, ref, aa=False):
    res = {}
    ambiguous = 'X' if aa else 'N'

    with open_file(aln_file, 'r') as fh:
        for si, (name, seq) in enumerate(SimpleFastaParser(fh)):
            if si%1000==0 and si:
                print(f"sequence {si}")
            res[name] = ",".join([f"{a}{p}{d}" for a,p,d in get_differences(ref, seq, ambiguous)])

    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="transform nextalign output to sparse format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--directory', type=str, required=True, help="directory with nextalign output")
    parser.add_argument('--alignment', type=str, required=False, help="nucleotide alignment (if not part of default pattern)")
    parser.add_argument('--insertions', type=str, required=False, help="insertions (if not part of default pattern)")
    parser.add_argument('--basename', type=str, required=True, help="output pattern")
    parser.add_argument('--reference', type=str, required=True, help="reference sequence")
    parser.add_argument('--genes', nargs="+", required=True, help="list of gene names to summarize mutations for")
    parser.add_argument('--genemap', type=str, required=True, help="annotation")
    parser.add_argument('--output', type=str, required=True, help="output tsv file")
    args = parser.parse_args()

    res = read_reference(args.reference, args.genemap)
    ref = res['nuc']
    translations = res['translations']

    nucleotide_alignment = args.alignment or os.path.join(args.directory, args.basename+'.aligned.fasta*')
    insertions = os.path.join(args.directory, args.basename+'.insertions.csv')

    genes = set(args.genes)
    gene_files = glob.glob(os.path.join(args.directory, args.basename+'.gene.*.fasta*'))

    compressed = {}
    res = to_mutations(nucleotide_alignment, ref)
    compressed = {'nucleotide_mutations': pd.DataFrame(res.values(), index=res.keys(), columns=['nucleotide'])}
    for fname in gene_files:
        # Find the gene name in the current gene file, since the filename may have multiple suffixes.
        gene = (set(fname.split('.')) & genes).pop()
        res = to_mutations(fname, translations[gene], aa=True)
        compressed[gene] = pd.DataFrame(res.values(), index=res.keys(), columns=[gene])

    res = pd.concat(compressed.values(), axis=1)
    res.to_csv(args.output, sep='\t')
