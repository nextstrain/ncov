"""
Mask initial bases from alignment FASTA
"""
import argparse, gzip, sys
sys.path.insert(0,'.')
from collections import defaultdict
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from priorities import sequence_to_int_array
from augur.utils import read_metadata
from datetime import datetime

tmrca = datetime(2019, 12, 1).toordinal()

def expected_divergence(date, rate_per_day = 25/365):
    try:
        return (datetime.strptime(date, '%Y-%m-%d').toordinal() - tmrca)*rate_per_day
    except:
        return np.nan


def analyze_divergence(sequences, metadata, reference, mask_5p=0, mask_3p=0):
    int_ref = sequence_to_int_array(reference, fill_gaps=False)
    diagnostics = defaultdict(dict)
    fill_value = 110
    gap_value = 45
    ws = 50
    known_true_clusters = [(28880,28883)]
    known_true_cluster_array = np.ones_like(int_ref, dtype=int)
    for b,e in known_true_clusters:
        known_true_cluster_array[b:e]=0

    cluster_cut_off = 10
    with open(sequences) as fasta:
        for h,s in SimpleFastaParser(fasta):
            left_gaps = len(s) - len(s.lstrip('-'))
            right_gaps = len(s) - len(s.rstrip('-'))
            s = sequence_to_int_array(s, fill_value=fill_value, fill_gaps=False)
            # mask from both ends to avoid exclusion for problems at sites that will be masked anyway
            if mask_5p:
                s[:mask_5p] = fill_value
            if mask_3p:
                s[-mask_3p:] = fill_value

            # fill terminal gaps -- those will be filled anyway
            if left_gaps:
                s[:left_gaps] = fill_value
            if right_gaps:
                s[-right_gaps:] = fill_value

            # determine non-gap non-N mismatches
            snps = (int_ref!=s) & (s!=fill_value) & (s!=gap_value)
            # determine N positions
            filled = s==fill_value
            # determine gap positions (cast to int to detect start and ends)
            gaps = np.array(s==gap_value, dtype=int)
            gap_start = np.where(np.diff(gaps)==1)[0]
            gap_end = np.where(np.diff(gaps)==-1)[0]

            # determined mutation clusters by convolution with an array of ones => running window average
            clusters = np.array(np.convolve(snps*known_true_cluster_array, np.ones(ws), mode='same')>=cluster_cut_off, dtype=int)
            # determine start and end of clusters. extend by half window size on both ends.
            cluster_start = [0] if clusters[0] else []
            cluster_start.extend([max(0, x-ws//2) for x in np.where(np.diff(clusters)==1)[0]])
            cluster_end = [min(int_ref.shape[0], x+ws//2) for x in np.where(np.diff(clusters)==-1)[0]]
            if clusters[-1]:
                cluster_end.append(int_ref.shape[0])

            diagnostics[h] = {'snps':list(np.where(snps)[0]), 'gaps': list(zip(gap_start, gap_end)), 'gap_sum':np.sum(gaps),
                              'no_data':np.sum(filled) - mask_3p - mask_5p,
                              'clusters': [(b,e,np.sum(snps[b:e])) for b,e in zip(cluster_start, cluster_end)]}

    return diagnostics

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", type=str, required=True, help="FASTA file of alignment")
    parser.add_argument("--reference", type = str, required=True, help="reference sequence")
    parser.add_argument("--metadata", type = str, required=True, help="metadata")
    parser.add_argument("--mask-from-beginning", type = int, default=0, help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type = int, default=0, help="number of bases to mask from end")
    parser.add_argument("--output-diagnostics", type=str, required=True, help="Output of stats for every sequence")
    parser.add_argument("--output-flagged", type=str, required=True, help="Output of sequences flagged for exclusion with specific reasons")
    parser.add_argument("--output-exclusion-list", type=str, required=True, help="Output to-be-reviewed addition to exclude.txt")
    args = parser.parse_args()

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    ref  = SeqIO.read(args.reference, 'genbank').seq
    metadata, _ = read_metadata(args.metadata)

    diagnostics = analyze_divergence(args.alignment, metadata, ref,
                                     mask_5p=args.mask_from_beginning,
                                     mask_3p=args.mask_from_end)
    snp_cutoff = 25
    no_data_cutoff = 3000
    flagged_sequences = []
    # output diagnostics for each sequence, ordered by divergence
    with open(args.output_diagnostics, 'w') as diag:
        diag.write('\t'.join(['strain', 'divergence', 'excess divergence', '#Ns', '#gaps', 'clusters', 'gaps', 'all_snps', 'gap_list'])+'\n')
        for s, d in sorted(diagnostics.items(), key=lambda x:len(x[1]['snps']), reverse=True):
            expected_div = expected_divergence(metadata[s]['date']) if s in metadata else np.nan
            diag.write('\t'.join(map(str,[s, len(d['snps']), round(len(d['snps']) - expected_div,2),
                     d['no_data'], d['gap_sum'],
                     ','.join([f'{b}-{e}' for b,e,n in d['clusters']]),
                     ','.join([f'{b}-{e}' for b,e in d['gaps']]),
                     ','.join(map(str, d['snps'])),
                     ",".join([",".join([str(x) for x in range(b,e)]) for b,e in d["gaps"]])]))+'\n')


            msg = ""
            reasons = []
            if not np.isnan(expected_div) and np.abs(len(d['snps']) - expected_div) > snp_cutoff:
                msg += f"too high divergence {np.abs(len(d['snps']) - expected_div):1.2f}>{snp_cutoff};"
                reasons.append('divergence')
            if len(d['clusters']):
                msg += f"{len(d['clusters'])} SNP clusters with {','.join([str(x[2]) for x in d['clusters']])} SNPs each;"
                reasons.append('clustered mutations')
            if d['no_data']>no_data_cutoff:
                msg += f"too many Ns ({d['no_data']}>{no_data_cutoff})"
                reasons.append('too many ambigous sites')

            if msg:
                flagged_sequences.append([s, msg, tuple(reasons), metadata.get(s,{})])

    # write out file with sequences flagged for exclusion sorted by date
    to_exclude_by_reason = defaultdict(list)
    with open(args.output_flagged, 'w') as flag:
        flag.write(f'strain\tcollection_date\tsubmission_date\tflagging_reason\n')
        for s, msg, reasons, meta in sorted(flagged_sequences, key=lambda x:x[3].get('date_submitted', 'XX'), reverse=True):
            flag.write(f"{s}\t{metadata[s]['date'] if s in metadata else 'XXXX-XX-XX'}\t{metadata[s]['date_submitted'] if s in metadata else 'XXXX-XX-XX'}\t{msg}\n")
            to_exclude_by_reason[reasons].append(s)

    # write out file with sequences flagged for exclusion sorted by date
    with open(args.output_exclusion_list, 'w') as excl:
        for reason in to_exclude_by_reason:
            excl.write(f'\n# {"&".join(reason)}\n')
            excl.write('\n'.join(to_exclude_by_reason[reason])+'\n')

