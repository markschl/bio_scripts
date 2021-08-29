#!/usr/bin/env python2

import pysam


def fix(input, output='-', outformat=''):

    infile = pysam.AlignmentFile(input, 'rb')
    outfile = pysam.AlignmentFile(output, 'w' + outformat, header=infile.header)

    for read in infile:
        cigar = read.cigar
        if cigar:
            # remove insertion at start
            if cigar[0][0] == 2:
                read.reference_start += cigar[0][1]
                cigar = cigar[1:]
            # remove insertion at end
            if cigar[-1][0] == 2:
                cigar = cigar[:-1]

            read.cigar = cigar
        outfile.write(read)


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser('''
        Modifies SAM files produced by VSEARCH (https://github.com/torognes/vsearch) in order to be easier to handle.
        Deletions at the start and insertions at the end are removed from the CIGAR, and start positions adjusted.
    ''')
    p.add_argument('input', help='')
    p.add_argument('-o', '--output', default='-', help='')
    p.add_argument('-f', '--outformat', default='', help='')
    args = p.parse_args()

    fix(**vars(args))
