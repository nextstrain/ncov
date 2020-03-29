"""
Mask initial bases from alignment FASTA
"""
import argparse
from random import shuffle
from collections import defaultdict
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

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    all_seqs_dict = {x.id: x.upper() for x in AlignIO.read(args.alignment, format = 'fasta')}
    focal_seqs_dict = {x.id:x.upper() for x in AlignIO.read(args.focal_alignment, format='fasta')}
    print("Done reading the alignments.")

    # convert focal sequences to a numpy masked array
    focal_sequences = list(focal_seqs_dict.keys())
    focal_seqs_array = np.ma.array([focal_seqs_dict[x] for x in focal_sequences])
    focal_seqs_array.mask = np.bitwise_or(focal_seqs_array.data=='N', focal_seqs_array.data=='-')
    fraction_masked_in_focal = focal_seqs_array.mask.mean(axis=1)

    # remove focal sequences from all sequence to provide context data set
    context_seqs = [v for k,v in all_seqs_dict.items() if k not in focal_seqs_dict]
    # convert context set to masked arry
    context_seq_array = np.ma.array(context_seqs)
    context_seq_array.mask = np.bitwise_or(context_seq_array.data=='N', context_seq_array.data=='-')
    print("Done converting the alignments.")
    # for each context sequence, calculate minimal distance to focal set
    minimal_distance_to_focal_set = {}
    for i, seq in enumerate(context_seqs):
        closest_match = (np.ma.sum(focal_seqs_array!=context_seq_array[i], axis=1) 
                         + fraction_masked_in_focal).argmin()
        min_dist = np.ma.sum(focal_seqs_array[closest_match]!=context_seq_array[i])
        if i%100==0:
            print(f"compared {i} out of {len(context_seqs)} sequences")
        # for each context sequence, store distance and index of closest focal match
        minimal_distance_to_focal_set[seq.id] = (min_dist, closest_match)

    # for each focal sequence with close matches (using the index), we list all close contexts
    close_matches = defaultdict(list)
    for seq in minimal_distance_to_focal_set:
        close_matches[minimal_distance_to_focal_set[seq][1]].append(seq)

    for f in close_matches:
        shuffle(close_matches[f])
        close_matches[f].sort(key=lambda x:minimal_distance_to_focal_set[x][0])

    # export priorities
    with open(args.output, 'w') as fh:
        for i, seq in enumerate(context_seqs):
            # use distance as negative priority
            # penalize masked (N or -) -- 333 masked sites==one mutations
            # penalize if many sequences are close to the same focal one by using the index of the shuffled list of neighbours
            # currently each position in this lists reduced priority by 0.2, i.e. 5 other sequences == one mutation
            position = close_matches[minimal_distance_to_focal_set[seq.id][1]].index(seq.id)
            priority = -minimal_distance_to_focal_set[seq.id][0] - 0.003*np.sum(context_seq_array[i].mask) - 0.2*position
            fh.write(f"{seq.id}\t{priority:1.2f}\n")
