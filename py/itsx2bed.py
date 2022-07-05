#!/usr/bin/env python3

import argparse
from collections import OrderedDict
import os
import re
import csv
from sys import stderr


class DomainWriter(object):
    def __init__(self, name, out_prefix, complementary_info=None):
        filename = '{}.{}.bed'.format(out_prefix, name.replace('.', '_'))
        f = open(filename, 'w')
        self.writer = csv.writer(f, delimiter='\t', lineterminator=os.linesep)
        self.n = 0
        self.name = name
        self.complementary_info = complementary_info
    
    def write(self, seq_id, start, end, length):
        self.n += 1
        row = [seq_id, start, end]
        if self.complementary_info is not None:
            if self.complementary_info[seq_id]:
                strand = '-'
                if reverse:
                    # reverse strand: coordinates in positions file are reversed as well.
                    # However, the query sequences are not reversed, so we have to adjust.
                    start, end = length - end, length - start
                    row[1] = start
                    row[2] = end
            else:
                strand = '+'
            row += [seq_id, '0', strand]
        if not inverted:
            self.writer.writerow(row)
        else:
            if start > 0:
                row[1] = 0
                row[2] = start
                self.writer.writerow(row)
            if end < length:
                row[1] = end
                row[2] = length
                self.writer.writerow(row)
            self.writer.writerow(row)


def parse_domain(d):
    m = domain_p.search(d)
    rng = m.group(2)
    domain = m.group(1)
    if rng in ('Not found', 'No end', 'No start'):
        return domain, None
    else:
        start, end = rng.split('-', 1)
    return domain, (int(start), int(end))



p = argparse.ArgumentParser(description='''
Reads the output of ITSx (http://microbiology.se/software/itsx/) and creates BED files for each of
the domains.
''')
p.add_argument('prefix')
p.add_argument('-o', '--out-prefix', help='Output prefix. Default: same as input prefix')
p.add_argument('-d', '--domains', default='SSU,ITS1,5.8S,ITS2,LSU,ITS',
               help='Comma delimited list of domains, for which BED files should be produced, and which should be'
                    'included in --domain-list')
p.add_argument('-l', '--domain-list', action='store_true',
               help='Output a tab-separated file containing a comma delimited list of domains for each sequence ID'
                    '(to <out_prefix>.domains.txt)')
p.add_argument('--minlen', type=int, default=4, 
               help='Minimum domain length for BED coordinates to be written')
p.add_argument('-p', '--partial-58S', action='store_true',
               help='Should partial 5.8S domains be included even though ITSx reports that the start/end positions'
                    'were not found? The domains will be included if the remaining sequence is at most --maxlen-58S'
                    'bp long. Otherwise, it is assumed to be broken (see --broken-58S). It usually makes sense to activate this.')
p.add_argument('-b', '--broken-58S', action='store_true',
               help='If 5.8S start/end was not found in the seuqence and the putative 5.8S seems longer than --maxlen-58S,'
                    'should an approximate 5.8S domain still be returned? It will be exactly --maxlen-58S bp wide, and all'
                    'remaining sequence is ignored.')
p.add_argument('--len-58S', type=int, default=200, 
               help='Maximum length that 5.8S sequences can have if inferred indirectly (in --approx-58S mode).'
                    '[default: 200]. Note that if the sequence portion after the ITS1 end / before the ITS2 start'
                    'is longer than --len-58S, it is possible that some ITS1 / ITS2 part is dropped completely.')
p.add_argument('-e', '--exclude', action='append', default=[], 
               help='Exclude if given string is found in comment of positions file.')
p.add_argument('--required-domains', 
               help='Comma delimited list of domains that are required to be found. Otherwise, the sequence will be excluded.')
p.add_argument('-r', '--reverse', action='store_true',
               help='Reverse the coordinate positions relative to the sequence if on the reverse strand.')
p.add_argument('-i', '--inverted', action='store_true',
               help='Additionally write inverted BED output: instead of returing the domain coordinates, output the inverse,'
                    'which are the regions before and after the domain *not* covered by the domain.')
p.add_argument('-v', '--verbose', action='store_true', help='Output some debug information.')
args = p.parse_args()


domain_p = re.compile(r'([^:]+): (.+)')
len_p = re.compile(r'(\d+) bp\.')
exclude = set(args.exclude)
required_domains = set(d.strip() for d in args.required_domains.split(',')) if args.required_domains is not None else None

minlen = args.minlen
prefix = args.prefix
s58_len = args.len_58S
s58_partial = args.partial_58S
s58_broken = args.broken_58S
inverted = args.inverted
verbose = args.verbose
reverse = args.reverse
out_prefix = prefix if args.out_prefix is None else args.out_prefix
posfile = prefix + '.positions.txt'
details = prefix + '.extraction.results'
domain_list = csv.writer(open(out_prefix + '.domains.txt', 'w'), delimiter='\t') if args.domain_list else None
out_domains = [d.strip() for d in args.domains.split(',')]

# Check if details file present
if os.path.exists(details):
    with open(details) as f:
        complementary_info = {l[0]: l[3] == '1' for l in csv.reader(f, delimiter='\t') if l}
else:
    complementary_info = None
    print('No details file found. Re-run with --detailed_results if there can be'
          'hits on the reverse strand', file=stderr)


# domains in order found in ITSx positions file, with 'ITS' as last -> order of writers
domains = ['SSU', 'ITS1', '5.8S', 'ITS2', 'LSU'] + ['ITS']

domain_writers = [None]*6

for name in out_domains:
    i = domains.index(name)
    domain_writers[i] = DomainWriter(name, prefix, complementary_info)

# Total number of lines
n = 0
partial_58s = 0
strange_58s = 0

# Errors that lead to a line being excluded
errors = ('ITS region too long!', 'Chimeric!')

with open(posfile) as f:
    for line in csv.reader(f, delimiter='\t'):
        n += 1
        comment = line[7]
        # ignore errors / excluded comments
        if not line or any(e in comment for e in errors) or any(e in comment for e in exclude):
            if verbose:
                print('Sequence excluded:\n{}', '\t'.join(line), file=stderr)
            continue
        # parse line
        length = int(len_p.search(line[1]).group(1))
        seq_id = line[0]
        pos = 1
        data = [parse_domain(d) for d in line[2:7]]
        # 5.8S start / end can be missing, then we infer them 
        # from ITS1 end / ITS2 start
        if data[2][1] is None and 'only partial 5.8S' in comment:
            partial_58s += 1
            if s58_partial:
                its1_end = None if data[1][1] is None else data[1][1][1]
                its2_start = None if data[3][1] is None else data[3][1][0]
                strange = False
                if its1_end is not None:
                    start = its1_end + 1
                elif its2_start is not None:
                    start = max(1, its2_start - s58_len)
                    if start > 1: 
                        strange = True
                if its2_start is not None:
                    end = its2_start
                else:
                    end = min(length, start + s58_len + 1)
                    if end < length: 
                        strange = True
                data[2] = (data[2][0], (start, end))
                if strange:
                    strange_58s += 1
                    if not s58_broken:
                        data[2] = (domain, None)
                if verbose:
                    print('5.8S domain coordinates inferred ({})\n  "{}" (length: {})\n  -> {}'.format(
                        comment,
                        '  '.join(line[3:6]), length,
                        '  '.join('{}: {}'.format(n, rng) for n, rng in data[1:4])
                        ), file=stderr)
                    if strange:
                        print('  ... Parts of the sequence are not included!', file=stderr)
        
        # check if seq. should be excluded because of missing domains
        if required_domains is not None and \
            any(domain in required_domains and rng is None for domain, rng in data):
            if verbose:
                print('Not all domains were found:\n{}', '\t'.join(line), file=stderr)
            continue
        
        # write all domains to files
        global_start = None
        global_end = None
        for i in range(5):
            domain, rng = data[i]
            if rng is None:
                continue
            
            (start, end) = rng

            if i <= 3:
                if global_start is None and i >= 1:
                    global_start = start
                global_end = end

            # check length
            if end - start < minlen:
                continue
            # convert to 0-based start
            start -= 1

            # write to file if writer present
            if domain_writers[i] is not None:
                domain_writers[i].write(seq_id, start, end, length)

        # write BED for complete ITS
        if global_start is not None:
            assert global_end is not None
            if domain_writers[5] is not None:
                domain_writers[5].write(seq_id, global_start, global_end, length)

        # write domain list
        if domain_list is not None:
            detected = (d[0] for i, d in enumerate(data) if d[1] is not None and domain_writers[i] is not None)
            domain_list.writerow((seq_id, ','.join(detected)))


# write stats
print('Sequences from positions file written to BED\n(those not found at all are not included in the total count)', file=stderr)
for w in domain_writers:
    if w is not None:
        print('{}: {} of {} ({:.1f}%)'.format(w.name, w.n, n, w.n / n * 100), file=stderr)
if partial_58s > 0:
    msg = '' if s58_partial else ' (not included, use -p/--partial-58S to include)'
    print('* {} of {} sequences ({:.1f}%) had a partial 5.8S sequence{}'.format(partial_58s, n, partial_58s / n * 100, msg))
if strange_58s > 0:
    msg = '' if s58_broken else ' (not included, use -b/--broken-58S to include)'
    print('* {} of {} sequences ({:.1f}%) may have a broken 5.8S (start or end longer than {} bp, see --len-58S))'.format(strange_58s, n, strange_58s / n * 100, s58_len))
