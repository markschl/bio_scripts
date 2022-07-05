#!/usr/bin/env python3

import csv
import argparse
import re


def extract(taxonomy, output, ranks='', illegal=' ', unknown=None):
    ranks = ranks.split(',')
    writer = csv.writer(output, delimiter='\t')
    writer.writerow(['accession'] + ranks)
    illegal_patt = re.compile(r'[{}]'.format(re.escape(illegal)))
    if unknown is not None:
        unknown = re.compile(r'({})'.format('|'.join(re.escape(u) for u in unknown.split(','))), re.I)
    order = RankOrder()

    ranks = list(reversed(ranks))
    for row in csv.reader(taxonomy, delimiter='\t'):
        if not row:
            continue
        id = row[0]
        lineage = [n.split(':', 1) for n in reversed(row[1:])]
        lineage_dict = {r[0]: (r[1], i) for i, r in enumerate(lineage)}
        l = len(lineage)
        out = []
        # ranks are added in reverse order (from lowest to highest)
        # position within input lineage
        i = 0
        # position of next taxon that has no 'unknown' keyword
        known_i = None
        for r, rank in enumerate(ranks):
            if rank not in lineage_dict:
                # rank not present -> try to infer the best parent rank based on the rank name
                for j, x in enumerate(lineage[i+1:]):
                    rank2, name2 = x
                    cmp = order.cmp_names(rank, rank2)
                    if cmp is not None:
                        # matching rank found -> set position depending on relationship
                        i += j + cmp
                if i + 1 >= l:
                    name = 'unknown'
                else:
                    rank2, name2 = lineage[i + 1]
                    name = '{}_{}'.format(name2, rank2)
            else:
                name, i = lineage_dict[rank]

            if unknown is not None:
                # in case there are unknown names, proceed to higher ranks
                if known_i is not None and i >= known_i:
                    # all ranks up to this one were unknown, this one is known
                    known_i = None
                elif known_i is not None:
                    # current rank unknown (known_i is higher than i)
                    name = '{}_unknown'.format(lineage[known_i][1])
                else:
                    # check if unknown
                    j = i
                    while unknown.search(name) is not None:
                        j += 1
                        known_i = j
                        if j >= l:
                            name = 'unknown'
                            known_i = None
                            break
                        name = lineage[known_i][1]

                    if known_i is not None:
                        name = '{}_unknown'.format(name)

            if illegal:
               name = illegal_patt.sub('_', name)
            out.append(name)
        out.append(id)
        writer.writerow(reversed(out))


# species group / subgroup are not recognized separately
class RankOrder(object):
    """
    Compares taxonomic ranks
    """
    sub_pattern = re.compile(r'(sub|infra|parv)(?P<name>.+)')
    super_pattern = re.compile(r'super(?P<name>.+)')
    group_pattern = re.compile(r'(?P<name>[^ ]+) (sub)?group')

    def __init__(self):
        self.searched = {}

    def cmp_names(self, name1, name2):
        """
        @returns: None if the ranks are not related, 0 if they are equal, -1 if the second
            name is lower than the first, 1 if it is higher
        """
        #print('   cmp ', name1, name2)
        cmp = self.searched.get((name1, name2), 0)
        if cmp != 0:
            return cmp

        i1, i2 = 0, 0
        s = self.super_pattern.search(name1) or self.group_pattern.search(name1)
        if s is not None:
            name1 = s.group('name')
            i1 -= 1
        s = self.super_pattern.search(name2) or self.group_pattern.search(name2)
        if s is not None:
            name2 = s.group('name')
            i1 -= 1
        s = self.sub_pattern.search(name1)
        if s is not None:
            name1 = s.group('name')
            i1 += 1
        s = self.sub_pattern.search(name2)
        if s is not None:
            name2 = s.group('name')
            i2 += 1
        if name1 == name2:
            if i1 > i2:
                cmp = 1
            elif i1 == i2:
                cmp = 0
            else:
                cmp = -1
        else:
            cmp = None
        self.searched[(name1, name2)] = cmp
        return cmp


p = argparse.ArgumentParser("""
This script takes the output of 'get_taxonomy.py' and turns it into a tab delimited file with
names for the given ranks. If they are not found in the taxonomy, the script will try to find
the best name to return (based on sub/super hierarchies if present) and append its rank name.
This is done since there are many 'no_rank' classifications that would otherwise get lost,
and ranks are not always comparable between kingdoms.
Names containing one of the '--unknown' keywords will be replaced by the next defined name
upwards in the hierarchy, and '_unknown' is appended.
Make sure to always specify the ranks in the correct order, otherwise there will be confusion.
""")
p.add_argument("taxonomy", type=argparse.FileType('r'), help="File with Genbank accession on each line")
p.add_argument("-r", "--ranks", help="comma delimited list of ranks to extract",
               default="superkingdom,kingdom,phylum,class,order,family,genus,species,varietas")
p.add_argument("-o", "--output", type=argparse.FileType('w'), default='-')
p.add_argument("-u", "--unknown", default='environmental,uncultured,incertae,unidentified,unknown,unclassified')
p.add_argument("--illegal", default=' ', help='Illegal characters to replace by _')
args = p.parse_args()

extract(**vars(args))
