"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True,
                        help="FASTA file of alignment")
    parser.add_argument("--mask-from-beginning", type=int, required=True,
                        dest='begin_length', default=0,
                        help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type=int,
                        dest='end_length', default=0,
                        help="number of bases to mask from end")
    parser.add_argument("--mask-sites", nargs='+', type=int, default=[],
                        help="list of (1-based) sites to mask")
    parser.add_argument("--output", required=True,
                        help="FASTA file of output alignment")
    args = parser.parse_args()

    start = ["N"] * args.begin_length
    end = ["N"] * args.end_length

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            middle = list(seq[args.begin_length:-args.end_length])
            seq_list = start + middle + end
            for site in args.mask_sites:
                seq_list[site - 1] = "N"
            record.seq = Seq("".join(seq_list))
            Bio.SeqIO.write(record, outfile, 'fasta')
