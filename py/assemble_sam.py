#!/bin/env python3

import os
import sys
import re
import subprocess
from tempfile import NamedTemporaryFile
from io import StringIO
import multiprocessing

import pysam

from lib import FastaIO, seq
from lib.seq import SeqRecord
from lib.alignment import consensus


MAFFT = 'mafft'
FNULL = open(os.devnull, 'w')


def assemble(sam, reference, outfile=sys.stdout, key=None, processes=1, align_prefix=None, output_empty=False):
    if key is not None:
        keyobj = [Key(k) for k in key]
    else:
        keyobj = []

    samfile = pysam.Samfile(sam)

    pool = multiprocessing.Pool(processes)

    refseqs = {}
    refs = set(samfile.references)
    for rec in FastaIO.from_file(reference):
        if rec.id in refs:
            refseqs[rec.id] = rec
            refs.remove(rec.id)
    if refs:
        print("Not all references were found; {}".format(', '.join(refs)), file=sys.stderr)

    results = []

    for ref in samfile.references:
        if ref not in refseqs:
            continue

        keymap = {}
        for segment in samfile.fetch(ref):
            key = ''.join(k.get_key(segment.query_name) for k in keyobj)
            if key not in keymap:
                seqmap = keymap[key] = {}
            else:
                seqmap = keymap[key]
            seq = segment.query_sequence
            if seq not in seqmap:
                seqmap[seq] = [1, segment.get_tag('NM'), segment.query_name, key]
            else:
                dat = seqmap[seq]
                dat[0] += 1
                if segment.get_tag('NM') != dat[1]:
                    print('not same', segment.query_name, file=sys.stderr)

        refseq = refseqs[ref].seq
        if keymap:
            #print >> sys.stderr, ref
            chosen = [(refseq, (1, 0, ref, 'ref'))]
            for key, seqmap in keymap.items():
                # sort by (key, edit_distance, -size)  [biggest size first, then smallest edit distance]
                seqs = sorted(seqmap.items(), key=lambda x: (0 - x[1][0], x[1][1]))
                chosen.append(seqs[0])
            res = pool.apply_async(get_consensus, (ref, chosen, align_prefix))
            results.append(res)
            #FastaIO.write(get_consensus(ref, chosen, align_prefix), outfile)
        elif output_empty:
            FastaIO.write(SeqRecord(ref, refseq), outfile)

    for res in results:
        res.wait()
        FastaIO.write(res.get(), outfile)


_invalid_chars = re.compile('[{}]'.format(re.escape('<>:"/\\|?*')))


def get_consensus(ref, seqdata, align_prefix):

    seqs = [('{};{}x;dist={}'.format(key, count, distance), seq)
            for (seq, (count, distance, seqid, key)) in seqdata]

    f = NamedTemporaryFile(mode='wt')
    FastaIO.write((SeqRecord(id, seq) for id, seq in seqs), f)
    f.flush()

    aln = str(subprocess.check_output([MAFFT, f.name], stderr=FNULL), 'utf-8')
    cons = consensus(FastaIO.parse(StringIO(aln)))

    if align_prefix:
        FastaIO.to_file(FastaIO.parse(StringIO(aln)), '{}{}.fasta'.format(align_prefix, _invalid_chars.sub('_', ref)))
        f.close()

    cons = cons.replace('-', '')
    return SeqRecord('{} {}'.format(ref, ' '.join(s[0] for s in seqs)), cons)


# knows how to search key=value in a string
class Key(object):
    def __init__(self, key):
        self.patt = re.compile(";%s=([^;]+)" % re.escape(key))

    def get_key(self, string):
        m = self.patt.search(string)
        if m is None:
            return ''
        return m.group(1)


if __name__ == '__main__':
    from argparse import ArgumentParser, FileType

    p = ArgumentParser('''
    Extracts sequences from a BAM/SAM alignment and assembles them to one.
    The most abundant sequences from each amplicon group will be chosen
    (assumed to be the 'true' OTU sequences) and aligned using MAFFT.
    In case of ties, the sequence with the smallest edit distance to the OTU
    sequence is selected.
    The group IDs are expected to be found in the sequence ID, e.g.:
    'some_seq123;group=A'. The name of the group key can be specified using -k
    ('-k group' in this example).
    ''')
    p.add_argument("sam", help='Path to SAM/BAM file. The type is auto-recognized.')
    p.add_argument("reference", help='Path to SAM reference sequences')
    p.add_argument("-o", "--outfile", default=sys.stdout, type=FileType('w'),
                   help="Output file name or '-' for standard output (stdout). Defaults to stdout")
    p.add_argument("-d", "--align_prefix",
                   help='If specified, alignemnt files will be written to a location with the given prefix (with their id appended)')
    p.add_argument("-s", "--output_empty", action='store_true',
                   help='Specify if empty references found in the file should be written to output')
    p.add_argument("-k", "--key", action='append',
                   help="Key that specifies the name of the primer combination")
    p.add_argument("-p", "--processes", type=int, default=1,
                   help='Number of separate processes to use')
    args = p.parse_args()

    assemble(**vars(args))
