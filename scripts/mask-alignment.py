"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq

def mask_terminal_gaps(seq):
    L = len(seq)
    seq_trimmed = seq.lstrip('-')
    left_gaps = L - len(seq_trimmed)
    seq_trimmed = seq_trimmed.rstrip('-')
    right_gaps = L - len(seq_trimmed) - left_gaps
    return "N"*left_gaps + seq_trimmed + "N"*right_gaps 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument("--mask-terminal-gaps", action='store_true', help="fill all terminal gaps with N as they likely represent missing data")
    parser.add_argument("--mask-from-beginning", type = int, required=True, help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type = int, help="number of bases to mask from end")
    parser.add_argument("--mask-sites", nargs='+', type = int,  help="list of sites to mask")
    parser.add_argument("--output", required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    begin_length = 0
    if args.mask_from_beginning:
        begin_length = args.mask_from_beginning
    end_length = 0
    if args.mask_from_end:
        end_length = args.mask_from_end

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            if args.mask_terminal_gaps:
                seq = mask_terminal_gaps(seq)

            start = "N" * begin_length
            middle = seq[begin_length:-end_length]
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if args.mask_sites:
                for site in args.mask_sites:
                    seq_list[site-1] = "N"
            record.seq = Seq("".join(seq_list))
            Bio.SeqIO.write(record, outfile, 'fasta')
