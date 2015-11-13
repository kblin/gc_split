#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

_version = '0.1.0'


def main():
    """parse the command line and run"""
    parser = argparse.ArgumentParser(description="Split a set of FASTA records into high-GC and low-GC records")

    parser.add_argument('infiles', metavar='FILE', nargs='+',
                        help='File containing one or more sequence records')
    parser.add_argument('--high', default='high_gc.fa',
                        help='Filename for high GC records')
    parser.add_argument('--low', default='low_gc.fa',
                        help='Filename for low GC records')
    parser.add_argument('-t', '--threshold', type=float, default=50.0,
                        help='GC percentage over which to consider records as high GC (default: 50)')
    parser.add_argument('--format', default='fasta',
                        help='Select what fileformat the input files are in')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(_version))

    args = parser.parse_args()

    records = load_records(args.infiles, args.format)
    high_gc, low_gc = split_by_gc(records, args.threshold)
    num_high = SeqIO.write(high_gc, args.high, args.format)
    num_low = SeqIO.write(low_gc, args.low, args.format)

    print("Wrote {} high GC and {} low GC records".format(num_high, num_low))



def load_records(infiles, format):
    """Load records specified in infiles and concatenate them"""

    records = []

    for infile in infiles:
        records.extend(list(SeqIO.parse(infile, format)))

    return records


def split_by_gc(records, threshold=50):
    """Split records by DNA GC content.
    All seqs with GC > threshold will go into high_gc, others into low_gc
    """

    high_gc = []
    low_gc = []

    for record in records:
        if GC(record.seq) > threshold:
            high_gc.append(record)
        else:
            low_gc.append(record)

    return high_gc, low_gc


if __name__ == '__main__':
    main()
