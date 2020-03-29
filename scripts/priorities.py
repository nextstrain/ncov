"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import numpy as np
import Bio.SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", type=str, required=True, help="FASTA file of alignment")
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--focal-alignment", type = str, required=True, help="focal smaple of sequences")
    parser.add_argument("--output", type=str, required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    all_seqs_dict = {x.id: x.upper() for x in AlignIO.read(args.alignment, format = 'fasta')}
    focal_seqs_dict = {x.id:x.upper() for x in AlignIO.read(args.focal_alignment, format='fasta')}

    focal_seqs_array = np.ma.array([x for x in focal_seqs_dict.values()])
    focal_seqs_array.mask = np.bitwise_or(focal_seqs_array.data=='N', focal_seqs_array.data=='-')

    context_seqs = [v for k,v in all_seqs_dict.items() if k not in focal_seqs_dict]
    context_seq_array = np.ma.array(context_seqs)
    context_seq_array.mask = np.bitwise_or(context_seq_array.data=='N', context_seq_array.data=='-')

    minimal_distance_to_focal_set = {}
    for i, seq in enumerate(context_seqs):
        min_dist = np.ma.sum(focal_seqs_array!=context_seq_array[i], axis=1).min()
        print(seq.id, min_dist)
        minimal_distance_to_focal_set[seq.id] = min_dist


    with open(args.output, 'w') as fh:
        for i, seq in enumerate(context_seqs):
            priority = -minimal_distance_to_focal_set[seq.id] - 0.003*np.sum(context_seq_array[i].mask) 
            fh.write(f"{seq.id}\t{priority:1.2f}\n")
