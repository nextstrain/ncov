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
            attributes = {x.split('=')[0]:'='.join(x.split('=')[1:]) for x in entries[8].split(';')}
            if 'gene_name' in attributes:
                name = attributes['gene_name'].strip('"')
            else:
                name = None
            translation = Seq.translate(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start-1, end, strand=-1 if strand=='-' else 1)).extract(ref))
            translations[name] = str(translation)

    return {"nuc":ref, "translations":translations}

def summarise_differences(ref, query, isAA):
    """
    Summarise the differences between a provided reference and a query
    (both of which are numpy arrays with dtype int8)
    Returns a string of comma-seperated mutations
    """
    # in int8:   A = 65       T = 84      C = 67       G = 71      N = 78       - = 45      X = 88
    ambiguous = 88 if isAA else 78 # 'X' or 'N'
    # replace all leading and trailing gaps with the ambiguous character
    idx_not_gaps = np.where(query!=45)[0] # 45 is '-' (gap)
    if idx_not_gaps.size:
        query[0:idx_not_gaps[0]] = ambiguous
        query[idx_not_gaps[-1]+1:len(query)] = ambiguous
    else:
        # the query is nothing but gaps! We don't report any mutations here
        return ""
    # sometimes the query length is longer than the reference. In this case we preserve previous behavior
    # by discarding extra characters in the query
    if query.size>ref.size:
        query = query[0:ref.size]
    # find indicies where the query differs from the reference, and is not ambiguous
    changes = np.logical_and(ref!=query, query!=ambiguous).nonzero()[0]
    # save these as a comma-seperated list of <from><base><to>, where the base (position) is 1-based
    return ",".join([f"{chr(ref[idx])}{idx+1}{chr(query[idx])}" for idx in changes])

def to_numpy_array(input_string):
    return np.frombuffer(input_string.upper().encode('utf-8'), dtype=np.int8).copy()

def to_mutations(aln_file, ref, aa=False):
    res = {}
    ref_array = to_numpy_array(ref)
    with open_file(aln_file, 'r') as fh:
        for name, seq in SimpleFastaParser(fh):
            res[name] = summarise_differences(ref_array, to_numpy_array(seq), aa)
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
