"""Extract sequences from a given FASTA file that match the given list of sample names.
"""
import numpy as np
import argparse, sys
from Bio import AlignIO, SeqIO, Seq, SeqRecord
from Bio.AlignIO import MultipleSeqAlignment
from augur.translate import safe_translate

scoring_params = {"score_match":3, "score_mismatch":-1, "score_gapext":-1, "score_gapopen":-10}

def align_pairwise(seq1, seq2):
    bandw = int(np.abs((len(seq1)-len(seq2)))) + 100
    try:
        from seqanpy import align_overlap
        return align_overlap(seq1, seq2, **scoring_params, band=bandw)
    except ImportError:
        from Bio import pairwise2
        aln = pairwise2.align.globalms(seq1, seq2,
            scoring_params['score_match'], scoring_params['score_mismatch'],
            scoring_params['score_gapopen'], scoring_params['score_gapext'],
            penalize_end_gaps=False, one_alignment_only=True)[0]
        return aln[2], aln[0], aln[1]

def ref_align(ref, seq):
    score, refaln, seqaln = align_pairwise(ref, seq)

    ref_aln_array = np.array(list(refaln))
    seq_aln_array = np.array(list(seqaln))

    # stip gaps
    ungapped = ref_aln_array!='-'
    ref_aln_array_ungapped = ref_aln_array[ungapped]
    seq_aln_array_ungapped = seq_aln_array[ungapped]
    return ''.join(seq_aln_array_ungapped)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extract sample sequences by name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--reference", required=True, help="annotated genbank file")
    parser.add_argument("--output", required=True, help="FASTA file of extracted sample sequences")
    args = parser.parse_args()

    aln = SeqIO.parse(args.sequences, 'fasta')
    ref = SeqIO.read(args.reference, 'genbank')
    refstr = str(ref.seq).upper()

    alignment = []
    for seq in aln:
        print(seq.id)
        seq_aln = ref_align(refstr, str(seq.seq).upper())
        if seq_aln:
            if len(seq_aln)!=len(refstr):
                print(seq.name, seq_aln, refstr)
            else:
                seq.seq=Seq.Seq(seq_aln)
                alignment.append(seq)

    # output
    AlignIO.write(MultipleSeqAlignment(alignment), args.output, 'fasta')
