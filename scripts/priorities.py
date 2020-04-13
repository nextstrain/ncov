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
from Bio import AlignIO
from scipy import sparse


INITIALISATION_LENGTH = 1000000

# Function adapted from https://github.com/gtonkinhill/pairsnp-python
def calculate_snp_matrix(fastafile, consensus=None, zipped=False):
    # This function generate a sparse matrix where differences to the consensus are coded as integers.

    row = np.empty(INITIALISATION_LENGTH)
    col = np.empty(INITIALISATION_LENGTH, dtype=np.int64)
    val = np.empty(INITIALISATION_LENGTH, dtype=np.int8)

    r = 0
    n_snps = 0
    nseqs = 0
    seq_names = []
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
                consensus = np.frombuffer(s.lower().encode('utf-8'), dtype=np.int8).copy()
                consensus[(consensus!=97) & (consensus!=99) & (consensus!=103) & (consensus!=116)] = 97
            else:
                align_length = len(consensus)

            nseqs +=1
            seq_names.append(h)

            if(len(s)!=align_length):
                raise ValueError('Fasta file appears to have sequences of different lengths!')

            s = np.frombuffer(s.lower().encode('utf-8'), dtype=np.int8).copy()
            s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
            snps = consensus!=s
            right = n_snps + np.sum(snps)

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

    return {'snps': sparse_snps, 'consensus': consensus, 'names': seq_names}

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
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--focal-alignment", type = str, required=True, help="focal smaple of sequences")
    parser.add_argument("--output", type=str, required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    context_seqs_dict = calculate_snp_matrix(args.alignment)
    focal_seqs_dict = calculate_snp_matrix(args.focal_alignment, consensus = context_seqs_dict['consensus'])
    alignment_length = len(context_seqs_dict['consensus'])
    print("Done reading the alignments.")

    # read in focal alignment to figure out number N & gap
    focal_seqs = [x.upper() for x in AlignIO.read(args.focal_alignment, format='fasta')]
    focal_seqs_array = np.array([list(str(x.seq)) for x in focal_seqs])
    
    #focal_seqs2 = {x.id:x.upper() for x in AlignIO.read(args.focal_alignment, format='fasta')}
    #focal_seqs_array2 = np.array([list(str(focal_seqs[x].seq)) for x in focal_seqs])
    mask_count_focal = np.sum(np.logical_or(focal_seqs_array=='N',focal_seqs_array=='-'), axis=1)

    # remove focal sequences from all sequence to provide context data set
    keep = [(i, name) for i, name in enumerate(context_seqs_dict['names']) if name not in set(focal_seqs_dict['names'])]
    context_seqs_dict['snps'] = context_seqs_dict['snps'].tocsr()[[i for i, name in keep],:].tocsc()
    context_seqs_dict['names'] = [name for i, name in keep]

    tmp_context_seqs = {x.id: x.upper() for x in AlignIO.read(args.alignment, format='fasta')}
    context_seqs_array = np.array([list(str(tmp_context_seqs[x[1]].seq)) for x in keep])
    mask_count_context = {k[1]:m for k,m in zip(keep, np.sum(np.logical_or(context_seqs_array=='N',context_seqs_array=='-'), axis=1))}

    # for each context sequence, calculate minimal distance to focal set, weigh with number of N/- to pick best sequence
    d = np.array(calculate_distance_matrix(context_seqs_dict['snps'], focal_seqs_dict['snps'], consensus = context_seqs_dict['consensus']))
    closest_match = np.argmin(d+mask_count_focal/alignment_length, axis=1)
    #closest_match = np.argmin(d, axis=1)[:,0]
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
