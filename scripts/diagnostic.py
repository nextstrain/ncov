"""
Mask initial bases from alignment FASTA
"""
import argparse, gzip
from collections import defaultdict
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
from scripts.priorities import sequence_to_int_array
from augur.utils import read_metadata
from datetime import datetime

tmrca = datetime(2019, 11,15).toordinal()

def expected_divergence(date):
    try:
        return (datetime.strptime(date, '%Y-%m-%d').toordinal() - tmrca)*25/365
    except:
        return np.nan

def analyze_divergence(sequences, metadata, reference):
    int_ref = sequence_to_int_array(reference, fill_gaps=False)
    diagnostics = defaultdict(dict)
    fill_value = 110
    gap_value = 45
    ws = 100
    cluster_cut_off = 5
    with open(sequences) as fasta:
        for h,s in SimpleFastaParser(fasta):
            left_gaps = len(s) - len(s.lstrip('-'))
            right_gaps = len(s) - len(s.rstrip('-'))
            s = sequence_to_int_array(s, fill_value=fill_value, fill_gaps=False)
            if left_gaps:
                s[:left_gaps] = fill_value
            if right_gaps:
                s[-right_gaps:] = fill_value
            snps = (int_ref!=s) & (s!=fill_value) & (s!=gap_value)
            filled = s==fill_value

            gaps = np.array(s==gap_value, dtype=int)
            gap_start = np.where(np.diff(gaps)==1)[0]
            gap_end = np.where(np.diff(gaps)==-1)[0]

            clusters = np.array(np.convolve(snps, np.ones(ws), mode='same')>cluster_cut_off, dtype=int)
            cluster_start = np.where(np.diff(clusters)==1)[0]
            cluster_end = np.where(np.diff(clusters)==-1)[0]
            diagnostics[h] = {'snps':list(np.where(snps)[0]), 'gaps': list(zip(gap_start, gap_end)), 'gap_sum':np.sum(gaps),
                              'no_data':np.sum(filled),
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
    parser.add_argument("--output-diagnostics", type=str, required=True, help="Output of stats for every sequence")
    parser.add_argument("--output-flagged", type=str, required=True, help="Output of sequences flagged for exclusion")
    args = parser.parse_args()

    # load entire alignment and the alignment of focal sequences (upper case -- probably not necessary)
    ref  = SeqIO.read(args.reference, 'genbank').seq
    metadata, _ = read_metadata(args.metadata)

    diagnostics = analyze_divergence(args.alignment, metadata, ref)
    snp_cutoff = 20
    no_data_cutoff = 3000
    with open(args.output_diagnostics, 'w') as diag, open(args.output_flagged, 'w') as flag:
        diag.write('\t'.join(['strain', 'divergence', 'excess divergence', '#Ns', '#gaps', 'clusters', 'gaps', 'all_snps'])+'\n')
        for s, d in sorted(diagnostics.items(), key=lambda x:len(x[1]['snps'])):
            expected_div = expected_divergence(metadata[s]['date']) if s in metadata else np.nan
            diag.write('\t'.join(map(str,[s, len(d['snps']), round(len(d['snps']) - expected_div,2),
                     d['no_data'], d['gap_sum'],
                     ','.join([f'{b}-{e}' for b,e,n in d['clusters']]),
                     ','.join([f'{b}-{e}' for b,e in d['gaps']]),
                     ','.join(map(str, d['snps']))]))+'\n')


            msg = ""
            if not np.isnan(expected_div) and len(d['snps']) - expected_div > snp_cutoff:
                msg += f"too high divergence {len(d['snps']) - expected_div:1.2f}>{snp_cutoff};"
            if len(d['clusters']):
                msg += f"snps clusters;"
            if d['no_data']>no_data_cutoff:
                msg += f"too many Ns ({d['no_data']}>{no_data_cutoff})"
            if msg:
                flag.write(f'{s}\t{msg}\n')

