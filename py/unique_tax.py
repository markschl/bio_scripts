#!/usr/bin/env python3

import re
import csv
from collections import defaultdict

import argparse

from lib import FastaIO


def uniqify(fasta, output, conflict_out=None):
    # regex pattern for UTAX format (could be adapted for others)
    patt = re.compile(r'(.+?);tax=(.+?);?$')
    tax_patt = re.compile(r'([a-z]):([^,]+)')
    outer_fmt = '{};tax={};'
    rank_fmt = '{}:{}'
    sep = ','

    if conflict_out is not None:
        conflict_out = csv.writer(conflict_out, delimiter='\t')
    # dictionary containing the corresponding tree branch for each name (per rank)
    # {rank_name: {taxon_name: <tree_branch>, ...}}
    name_dict = defaultdict(dict)
    # taxonomic tree
    # {taxon_name: {lower_name: {...}, ...}}
    tree = {}
    # dict containing
    # {rank_name: {taxon_name: tree_branch_dict, ...}}
    changed = defaultdict(dict)
    for rec in FastaIO.parse(fasta):
        m = patt.search(rec.id)
        id = m.group(1)
        l = m.group(2)
        lineage = tax_patt.findall(l)
        branch = tree
        corrected = False

        parent_name = ''
        for i, r in enumerate(lineage):
            rank, name = r
            # check if already changed
            changed_rank = changed[rank]
            n = changed_rank.get((parent_name, name), None)
            if n is not None:
                if conflict_out is not None:
                    conflict_out.writerow([id, rank, name, n])
                name = n
                lineage[i] = (rank, name)
                corrected = True

            # get branch for given taxon
            if name not in branch:
                branch[name] = {}
            branch = branch[name]

            # ensure that the name is not duplicated by checking name_dict & storing branch there
            if name not in name_dict[rank]:
                name_dict[rank][name] = branch

            # not the same branch -> rename
            elif branch != name_dict[rank][name]:
                # name
                new_name = '{}_({})'.format(name, parent_name)
                changed_rank[(parent_name, name)] = new_name
                lineage[i] = (rank, new_name)
                corrected = True
                if conflict_out is not None:
                    conflict_out.writerow([id, rank, name, new_name])

            parent_name = name

        if corrected:
            rec.id = outer_fmt.format(id, sep.join(rank_fmt.format(rank, name) for rank, name in lineage))

        FastaIO.write(rec, output)


p = argparse.ArgumentParser("""
Searches for duplicate names within a UTAX formatted file that are on the same
taxonomic level. The first occurrence in the file will be retained, all others
will have _(name) appended where name is the name of the next higher rank.
This avoids problems with SINTAX classification.
""")
p.add_argument("fasta", type=argparse.FileType(), help="FASTA file with UTAX taxonomy annotations")
p.add_argument("-o", "--output", type=argparse.FileType('w'), default='-')
p.add_argument("-c", "--conflict-out", type=argparse.FileType('w'))
args = p.parse_args()

uniqify(**vars(args))
