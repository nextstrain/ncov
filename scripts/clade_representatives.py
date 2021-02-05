import argparse, os, glob
from Bio import SeqIO, SeqFeature, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
from collections import defaultdict
import json

def get_common_mutations(strains, mutation_dict):
    cmuts = defaultdict(lambda: defaultdict(int))
    for strain in strains:
        for gene, muts in mutation_dict.get(strain, {}).items():
            if muts:
                for mut in muts.split(','):
                    m = mut[0], int(mut[1:-1]), mut[-1]
                    cmuts[gene][m]+=1

    consensus = {}
    for gene in cmuts:
        consensus[gene]=[]
        for pos in sorted(cmuts[gene], key=lambda x:x[1]):
            if cmuts[gene][pos]>0.7*len(strains):
                consensus[gene].append(pos)

    return consensus


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="pull out public early representatives of clades",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--mutations", type=str, required=True, help="mutation summary file")
    parser.add_argument("--index", type=str, required=True, help="sequence statistics file")
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--output", type = str, required=True, help="output")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')
    mutations = pd.read_csv(args.mutations, sep='\t', index_col=0).fillna('')
    index = pd.read_csv(args.index, sep='\t', index_col=0)
    Ns = {i:x.N for i,x in index.iterrows()}
    mutation_dict = {i: x.to_dict() for i,x in mutations.iterrows()}

    metadata_public = metadata[metadata.genbank_accession!='?']

    clades = metadata.Nextstrain_clade.unique()
    good_sequences = []


    reps = {}
    for clade in clades:
        reps[clade] = {}
        clade_subset = metadata_public[metadata_public.Nextstrain_clade==clade].sort_values(by='date')
        clade_subset = clade_subset[clade_subset.strain.apply(lambda x:Ns.get(x,3e4)<5e5)]
        if len(clade_subset):
            print(clade, [x for x in clade_subset[:5].strain], len(clade_subset))
            muts = get_common_mutations(clade_subset.strain, mutation_dict)
            reps[clade]["representatives"] = [(x.strain, x.genbank_accession) for i,x in clade_subset[:5].iterrows()]
            reps[clade]["mutations"] = muts
        else:
            print("no public high quality sequences for clade", clade)


    with open(args.output, 'w') as fh:
        json.dump(reps,fh)
