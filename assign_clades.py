#!/usr/bin/env python3
"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import numpy as np
import argparse, sys, os
from Bio import AlignIO, SeqIO, Seq, SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from augur.translate import safe_translate
from augur.clades import read_in_clade_definitions, is_node_in_clade
from augur.align import generate_alignment_cmd, AlignmentError
from augur.utils import run_shell_command, load_features


class tmpNode(object):
    def __init__(self):
        self.sequences = {}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign clades to sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of HA sequences")
    parser.add_argument("--output", type=str, default='clade_assignment.tsv', help="tsv file to write clade definitions to")
    parser.add_argument("--keep-temporary-files", action='store_true', help="don't clean up")
    parser.add_argument("--chunk-size", default=10, type=int, help="don't clean up")
    parser.add_argument("--nthreads", default=1, type=int, help="Number of threads to use in alignment")
    args = parser.parse_args()

    refname = f"config/reference.gb"
    features = load_features(refname)
    seqs = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(refname, 'genbank')
    clade_designations = read_in_clade_definitions(f"config/clades.tsv")

    log_fname = "clade_assignment.log"
    in_fname = "clade_assignment_tmp.fasta"
    out_fname = "clade_assignment_tmp_alignment.fasta"


    output = open(args.output, 'w')
    print('name\tclade\tparent clades', file=output)

    # break the sequences into chunks, align each to the reference, and assign clades one-by-one
    done = False
    while not done:
        # generate a chunk with chunk-size sequences
        chunk = []
        while len(chunk)<args.chunk_size and (not done):
            try:
                seq = seqs.__next__()
                chunk.append(seq)
            except StopIteration:
                done = True

        print(f"writing {len(chunk)} and the reference sequence to file '{in_fname}' for alignment.")
        with open(in_fname, 'wt') as fh:
            SeqIO.write(chunk + [ref], fh, 'fasta')

        # align using augur command structure
        cmd = generate_alignment_cmd('mafft', args.nthreads, None, in_fname, out_fname, log_fname)
        success = run_shell_command(cmd)

        if not success:
            raise AlignmentError("Error during alignment")

        alignment = AlignIO.read(out_fname, 'fasta')
        for seq in alignment:
            if seq.id==ref.id:
                continue
            # read sequence and all its annotated features
            seq_container = tmpNode()
            seq_str = str(seq.seq)
            seq_container.sequences['nuc'] = {i:c for i,c in enumerate(seq_str)}
            for fname, feat in features.items():
                if feat.type != 'source':
                    seq_container.sequences[fname] = {i:c for i,c in enumerate(safe_translate(feat.extract(seq_str)))}

            # for each clade, check whether it matches any of the clade definitions in the tsv
            matches = []
            for clade_name, clade_alleles in clade_designations.items():
                if is_node_in_clade(clade_alleles, seq_container, ref):
                    matches.append(clade_name)

            # print the last match as clade assignment and all others as ancestral clades
            # note that this assumes that clades in the tsv are ordered by order of appearence.
            # furthermore, this will only work if parent clades don't have definitions that exclude
            # child clades, i.e. positions can only be additive for this to work.
            if matches:
                print(f"{seq.description}\t{matches[-1]}\t{', '.join(matches[:-1])}", file=output)
            else:
                print(f"{seq.description}\t -- \t", file=output)

    output.close()

    if not args.keep_temporary_files:
        os.remove(log_fname)
        os.remove(in_fname)
        os.remove(out_fname)

