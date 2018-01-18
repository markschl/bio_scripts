#!/usr/bin/env python2

from __future__ import print_function

import csv
import pysam

from lib.seq import reverse_complement


def get_range(read, start, end):
    cigar = read.cigartuples
    while cigar[-1][0] == D:
        # deletions at the end are interferring with get_offset()
        del cigar[-1]
    read_start = read.reference_start
    start -= read_start
    end -= read_start
    o0 = get_offset(cigar, start)
    o1 = get_offset(cigar, end)
    return start + o0, end + o1


# CIGAR codes understood by the program
M=0
I=1
D=2
N=3

def get_offset(cigar, pos):
    ref_pos = 0
    offset = 0
    for op, n in cigar:
        if op == D:
            ref_pos += n
            offset -= n
        elif op == M or op == N:
            ref_pos  += n
        elif op == I:
            offset += n
        else:
            print('CIGAR operation not supported', op)
        if ref_pos >= pos:
            break
    return offset


def run(bam, bed, output, margin=0, log=None, refout=False, minlen=1):
    if log is not None:
        log = csv.writer(log, delimiter='\t')

    with pysam.AlignmentFile(bam, 'rb') as f:
        for line in csv.reader(bed, delimiter='\t'):
            if not line:
                continue
            id = line[0]
            start = int(line[1])
            end = int(line[2])
            try:
                strand = line[5]
            except IndexError:
                strand = '+'
            if strand == '+':
                reverse = False
            elif strand == '-':
                reverse = True
            else:
                print('Invalid strand specification: ', strand)
                exit()

            writer = Writer(minlen, margin)
            for i, read in enumerate(f.fetch(region=id)):
                if i == 0 and refout:
                    offset = read.reference_start
                    writer.write_trim(
                        read.reference_name,
                        # TODO: what do the inserted bytes signify? Doesn't work with Python 3 due to coding issues
                        read.get_reference_sequence().decode('ascii', 'ignore').encode('ascii').upper(),
                        offset + start, offset + end,
                        output,
                        reverse
                    )

                s, e = get_range(read, start, end)
                writer.write_trim(read.query_name, read.query_sequence, s, e, output, reverse)

            if log is not None:
                log.writerow((id, writer.n_written, writer.n))


class Writer(object):
    def __init__(self, minlen=1, margin=0):
        self.n = 0
        self.n_written = 0
        self.minlen = minlen
        self.margin = margin

    def write_trim(self, seqname, seq, start, end, out, reverse=False):
        self.n += 1
        maxend = len(seq) - self.margin
        s = max(self.margin, min(maxend, start))
        e = max(self.margin, min(maxend, end))

        if e - s >= self.minlen:
            seq = seq[s:e]
            if reverse:
                seq = reverse_complement(seq)
            self.n_written += 1
            if end > maxend:
                e = '{}]'.format(e)
            out.write('>{}:{}-{}\n{}\n'.format(seqname, s+1, e, seq))


if __name__ == '__main__':

    import argparse

    p = argparse.ArgumentParser('''
    Extracts all reads in a BAM file within the regions specified in the BED
    file to a FASTA file. They are trimmed to exactly the given intervals
    (however, the starts/ends of the read can be removed using -m
    ''')
    p.add_argument('bam', help='')
    p.add_argument('bed', type=argparse.FileType(), help='')
    p.add_argument('-o', '--output', type=argparse.FileType('w'), default='-')
    p.add_argument('-l', '--log', type=argparse.FileType('w'), default='-')
    p.add_argument('-r', '--refout', action='store_true',
        help="Output the trimmed reference sequence. CAUTION: Seems not to work \
        correctly at present. Requires the MD tag to be set."
    )
    p.add_argument('-m', '--margin', type=int, default=0,
        help='Remove n bases from the start and end of aligned reads if they \
        happen to be within the BED range')
    args = p.parse_args()

    run(**vars(args))
