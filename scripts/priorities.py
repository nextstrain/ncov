"""
Mask initial bases from alignment FASTA
"""
import argparse
from random import shuffle
from collections import defaultdict
import Bio
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from scipy import sparse


def compactify_sequences(sparse_matrix, sequence_names):
    sequence_groups = defaultdict(list)
    for s, snps in zip(sequence_names, sparse_matrix):
        ind = snps.nonzero()
        vals = np.array(snps[ind])
        if len(ind[1]):
            sequence_groups[tuple(zip(ind[1], vals[0]))].append(s)
        else:
            sequence_groups[tuple()].append(s)

    return sequence_groups

INITIALISATION_LENGTH = 1000000

def sequence_to_int_array(s, fill_value=110, fill_gaps=True):
    seq = np.frombuffer(str(s).lower().encode('utf-8'), dtype=np.int8).copy()
    if fill_gaps:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116)] = fill_value
    else:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116) & (seq!=45)] = fill_value
    return seq

# Function adapted from https://github.com/gtonkinhill/pairsnp-python
def calculate_snp_matrix(fastafile, consensus=None, zipped=False, fill_value=110):
    # This function generate a sparse matrix where differences to the consensus are coded as integers.

    row = np.empty(INITIALISATION_LENGTH)
    col = np.empty(INITIALISATION_LENGTH, dtype=np.int64)
    val = np.empty(INITIALISATION_LENGTH, dtype=np.int8)

    r = 0
    n_snps = 0
    nseqs = 0
    seq_names = []
    filled_positions = []
    current_length = INITIALISATION_LENGTH
    if zipped:
        fh = gzip.open(fastafile, 'rt')
    else:
        fh = open(fastafile, 'rt')
    with fh as fasta:
        for h,s in SimpleFastaParser(fasta):
            if consensus is None:
                align_length = len(s)
                # Take consensus as first sequence
                consensus = sequence_to_int_array(s, fill_value=fill_value)
            else:
                align_length = len(consensus)

            nseqs +=1
            seq_names.append(h)

            if(len(s)!=align_length):
                raise ValueError('Fasta file appears to have sequences of different lengths!')

            s = sequence_to_int_array(s, fill_value=fill_value)
            snps = (consensus!=s) & (s!=fill_value)
            right = n_snps + np.sum(snps)
            filled_positions.append(np.where(s==fill_value)[0])

            if right >= (current_length/2):
                current_length = current_length + INITIALISATION_LENGTH
                row.resize(current_length)
                col.resize(current_length)
                val.resize(current_length)

            row[n_snps:right] = r
            col[n_snps:right] = np.flatnonzero(snps)
            val[n_snps:right] = s[snps]
            r += 1
            n_snps = right
    fh.close()

    if nseqs==0:
        raise ValueError('No sequences found!')

    row = row[0:right]
    col = col[0:right]
    val = val[0:right]

    sparse_snps = sparse.csc_matrix((val, (row, col)), shape=(nseqs, align_length))

    return {'snps': sparse_snps, 'consensus': consensus, 'names': seq_names, 'filled_positions': filled_positions}

# Function adapted from https://github.com/gtonkinhill/pairsnp-python
def calculate_distance_matrix(sparse_matrix_A, sparse_matrix_B, consensus):

    n_seqs_A = sparse_matrix_A.shape[0]
    n_seqs_B = sparse_matrix_B.shape[0]

    d = (1*(sparse_matrix_A==97)) * (sparse_matrix_B.transpose()==97)
    d = d + (1*(sparse_matrix_A==99) * (sparse_matrix_B.transpose()==99))
    d = d + (1*(sparse_matrix_A==103) * (sparse_matrix_B.transpose()==103))
    d = d + (1*(sparse_matrix_A==116) * (sparse_matrix_B.transpose()==116))

    d = d.todense()

    n_comp = (1*(sparse_matrix_A==110) * ((sparse_matrix_B==110).transpose())).todense()
    d = d + n_comp

    temp_total = np.zeros((n_seqs_A, n_seqs_B))
    temp_total[:] = (1*(sparse_matrix_A>0)).sum(1)
    temp_total += (1*(sparse_matrix_B>0)).sum(1).transpose()

    total_differences_shared = (1*(sparse_matrix_A>0)) * (sparse_matrix_B.transpose()>0)

    n_total = np.zeros((n_seqs_A, n_seqs_B))
    n_sum = (1*(sparse_matrix_A==110)).sum(1)
    n_total[:] = n_sum
    n_total += (1*(sparse_matrix_B==110)).sum(1).transpose()

    diff_n = n_total - 2*n_comp
    d = temp_total - total_differences_shared.todense() - d - diff_n

    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="generate priorities files based on genetic proximity to focal sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", type=str, required=True, help="FASTA file of alignment")
    parser.add_argument("--reference", type = str, required=True, help="reference sequence")
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--focal-alignment", type = str, required=True, help="focal smaple of sequences")
    parser.add_argument("--output", type=str, required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    ref = sequence_to_int_array(SeqIO.read(args.reference, 'genbank').seq)
    context_seqs_dict = calculate_snp_matrix(args.alignment, consensus=ref)
    focal_seqs_dict = calculate_snp_matrix(args.focal_alignment, consensus = ref)
    alignment_length = len(ref)
    print("Done reading the alignments.")

    # calculate number of masked sites in either set
    mask_count_focal = np.array([len(x) for x in focal_seqs_dict['filled_positions']])
    mask_count_context = {s: len(x) for s,x in zip(context_seqs_dict['names'], context_seqs_dict['filled_positions'])}

    # for each context sequence, calculate minimal distance to focal set, weigh with number of N/- to pick best sequence
    d = np.array(calculate_distance_matrix(context_seqs_dict['snps'], focal_seqs_dict['snps'], consensus = context_seqs_dict['consensus']))
    closest_match = np.argmin(d+mask_count_focal/alignment_length, axis=1)
    print("Done finding closest matches.")

    minimal_distance_to_focal_set = {}
    for context_index, focal_index in enumerate(closest_match):
        minimal_distance_to_focal_set[context_seqs_dict['names'][context_index]] = (d[context_index, focal_index], focal_seqs_dict["names"][focal_index])

    # for each focal sequence with close matches (using the index), we list all close contexts
    close_matches = defaultdict(list)
    for seq in minimal_distance_to_focal_set:
        close_matches[minimal_distance_to_focal_set[seq][1]].append(seq)

    for f in close_matches:
        shuffle(close_matches[f])
        close_matches[f].sort(key=lambda x: minimal_distance_to_focal_set[x][0] + mask_count_context[x]/alignment_length)

    # export priorities
    with open(args.output, 'w') as fh:
        for i, seqid in enumerate(context_seqs_dict['names']):
            # use distance as negative priority
            # penalize masked (N or -) -- 333 masked sites==one mutations
            # penalize if many sequences are close to the same focal one by using the index of the shuffled list of neighbours
            # currently each position in this lists reduced priority by 0.2, i.e. 5 other sequences == one mutation
            position = close_matches[minimal_distance_to_focal_set[seqid][1]].index(seqid)
            priority = -minimal_distance_to_focal_set[seqid][0] - 0.1*position
            fh.write(f"{seqid}\t{priority:1.2f}\n")
