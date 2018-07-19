#!/usr/bin/env python3

import argparse
import os
import re
import csv


p = argparse.ArgumentParser(description='''
Reads the output of ITSx (http://microbiology.se/software/itsx/) and creates BED files for each of
the domains.
''')
p.add_argument("prefix")
p.add_argument("-d", "--domains", default='SSU,ITS1,5.8S,ITS2,LSU')
p.add_argument("--minlen", type=int, default=4, help='minimum length of sequence')
p.add_argument("-e", "--exclude", action='append', help='Exclude if given string is found in comment of positions file.')
args = p.parse_args()

domain_p = re.compile(r'([^:]+): (.+)')
len_p = re.compile(r'(\d+) bp\.')
exclude = set(args.exclude)

minlen = args.minlen
prefix = args.prefix
posfile = prefix + '.positions.txt'
details = prefix + '.extraction.results'

outfiles = {name: csv.writer(open(prefix + '.%s.bed' % name.replace('.', '_'), 'w'), delimiter='\t', lineterminator=os.linesep)
            for name in args.domains.split(',')}


domains = ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')
s58_maxlen = 190

errors = ('ITS region too long!', 'Chimeric!')

n = 0
written = [0] * 5

if os.path.exists(details):
    with open(details) as f:
        complementary = {l[0]: l[3] == '1' for l in csv.reader(f, delimiter='\t') if l}
else:
    complementary = None
    print('No details file found. Re-run with --detailed_results if there are hits on the reverse strand')

with open(posfile) as f:
    for line in csv.reader(f, delimiter='\t'):
        comment = line[7]
        if line and not any(e in comment for e in errors) and not any(e in comment for e in exclude):
            length = int(len_p.search(line[1]).group(1))
            id = line[0]
            pos = 1
            for i in range(5):
                m = domain_p.search(line[i + 2])
                name = m.group(1)
                text = m.group(2)
                if text != 'Not found' and name in outfiles:
                    if text in ('No start', 'No end'):
                        assert name == '5.8S' # always true?
                        # get end from ITS2
                        start, end = None, None
                        its1 = domain_p.search(line[3]).group(2)
                        its1_end = None if its1 == 'Not found' else int(its1.split('-', 1)[1])
                        its2 = domain_p.search(line[5]).group(2)
                        its2_start = None if its2 == 'Not found' else int(its2.split('-', 1)[0])
                        if its1_end is None:
                            if its2_start is not None and its2_start > s58_maxlen:
                                start = its2_start - s58_maxlen
                            else:
                                start = 1
                        else:
                            start = its1_end + 1
                        if its2_start is None:
                            end = min(length - start + 1, s58_maxlen)
                        else:
                            end = its2_start
                    else:
                        start, end = map(int, text.split('-', 1))

                    n += 1
                    pos = start

                    if end - start >= minlen:
                        if complementary is None:
                            outfiles[name].writerow([id, start - 1, end])
                        else:
                            strand = '-' if complementary[id] else '+'
                            outfiles[name].writerow([id, start - 1, end, '', '', strand])
                        written[i] += 1

for domain, n_written in zip(domains, written):
    if n_written:
        print('{}: {} of {} ({:.1f}%)'.format(domain, n_written, n, n_written / n * 100))
