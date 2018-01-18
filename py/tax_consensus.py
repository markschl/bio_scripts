#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from collections import defaultdict
from itertools import islice
import re

from lib import FastaIO
from lib.util import *



def sorted_freqs(l):
    freqs = defaultdict(int)
    for item in l:
        freqs[item] += 1
    return sorted(list(freqs.items()), key=lambda k: k[1], reverse=True)


def consensus_tax(
        taxonomy, uc,
        fasta=None, out_prefix='consens_tax',
        threshold=-2, freq_threshold=0.6, alternative_threshold=0.7,
        max_alternatives=3, exclude='',
        empty_kws=''
):

    # generate dict of clusters
    clusters = defaultdict(list)
    for row in csv.reader(uc, delimiter='\t'):
        type = row[0]
        if type == 'C':
            target = row[8]
            if target not in clusters:
                clusters[target] = []
        elif type == 'H':
            ident = float(row[3])
            query = row[8]
            target = row[9]
            clusters[target].append((query, ident))

    # load taxonomy into memory
    with open(taxonomy) as t:
        taxreader = csv.reader(t, delimiter='\t')
        header = next(taxreader)
        tax = {row[0]: row[1:] for row in taxreader}

    # writers
    if out_prefix is None:
        out_prefix = os.path.join(os.path.dirname(taxonomy), 'consensus_tax')

    good_writer = csv.writer(open(out_prefix + '_good.txt', 'w'), delimiter='\t')
    good_writer.writerow(header)

    strange_writer = csv.writer(open(out_prefix + '_strange.txt', 'w'), delimiter='\t')
    strange_writer.writerow(header)

    singleton_writer = csv.writer(open(out_prefix + '_singletons.txt', 'w'), delimiter='\t')
    singleton_writer.writerow(header)

    # prepare exclude pattern
    exclude_patt = re.compile(r'({})'.format('|'.join(re.escape(kw) for kw in exclude.split(','))), re.I)
    empty_patt = re.compile(r'({})'.format('|'.join(re.escape(kw) for kw in empty_kws.split(','))))

    dist_exp = 2
    ranks = header[1:]
    ntaxa = len(ranks)
    # weight_sum = sum(i**weight_exp for i in range(1, ntaxa + 1))
    # weights = [(ntaxa - i)**weight_exp/weight_sum for i in range(ntaxa)]

    # adjust threshold
    if threshold < 0:
        threshold = ntaxa + threshold

    # contains (category, lineage)
    categorized = {}

    # loop
    for id, cluster in clusters.items():

        # get taxonomy and exclude by keyword
        good_clusters = [
            (id2, ident, tax[id2])
            for id2, ident in cluster
            if id2 in tax
        ]
        # sort by identities
        cluster_tax = sorted(good_clusters, key=lambda k: k[1], reverse=True)

        if id in tax:
            centroid_lineage = tax[id]
        elif len(cluster_tax) > 0:
            centroid_lineage = cluster_tax[0][2]
        else:
            categorized[id] = ('not_found', None)
            continue

        cluster_tax = [(id2, ident, lineage) for id2, ident, lineage in cluster_tax
                       if all(exclude_patt.search(name) is None for name in lineage)]

        nseqs = len(cluster_tax)
        if nseqs == 0:
            # no duplicates or excluded
            singleton_writer.writerow([id] + centroid_lineage)
            if all(exclude_patt.search(name) is None for name in centroid_lineage):
                categorized[id] = ('singleton', centroid_lineage)
            else:
                categorized[id] = ('excluded', centroid_lineage)
            continue

        # get frequencies of names at each rank level
        freqs = [defaultdict(float) for _ in range(ntaxa)]
        for id2, ident, lineage in cluster_tax:
            for i, item in enumerate(zip(lineage, freqs)):
                name, f = item
                f[name] += (ident / 100) ** dist_exp
        # sort decreasing identities
        freqs = [sorted([(name, n, n / nseqs) for name, n in f.items()], key=lambda k: k[1], reverse=True)
                 for f in freqs]

        # get lowest rank number where centroid taxonomy name
        # does not match the most abundant name of the cluster
        rank = 0
        for rank, item in enumerate(zip(centroid_lineage, freqs)):
            name, f = item
            if name != f[0][0]:
                break


        # check if similarity is enough high
        # TODO: simple threshold is not very clever
        if rank + 1 >= threshold:
            good_writer.writerow([id] + centroid_lineage)
            categorized[id] = ('good', centroid_lineage)
            continue

        #print('...\n{}\tcentroid (rank_f: {}, conf: {})\n{}'.format('\t'.join(centroid_lineage), f, freqs[rank][0][1]/nseqs, '\t'.join(f[0][0] for f in freqs)))

        # problematic -> propose another lineage
        # add names until their frequency drops too much
        proposed_lineage = []
        proposed_lineage_simple = []
        i = 0
        t = nseqs * freq_threshold
        for i, f in enumerate(freqs):
            first = f[0]
            if first[1] < t:
                break
            proposed_lineage.append(first[0])
            proposed_lineage_simple.append(first[0])

        for freq in freqs[i+1:]:
            first = freq[0]
            names = [first[0]]
            fsum = first[1]
            for j, item in islice(enumerate(freq), 1, max_alternatives):
                name, f, rel_freq = item
                fsum += f
                if fsum >= t:
                    break
                names.append(name)
            proposed_lineage.append('/'.join(names))
            proposed_lineage_simple.append(names[0])

        strange_writer.writerow([id] + centroid_lineage + [len(cluster), rank, freqs[rank][0][1]/nseqs])
        strange_writer.writerow([id + '_fix'] + proposed_lineage)
        categorized[id] = ('strange', proposed_lineage_simple)

    if fasta is not None:
        good_fasta = out_prefix + '_good.fasta'
        strange_fasta = out_prefix + '_strange.fasta'
        singleton_fasta = out_prefix + '_singletons.fasta'
        excluded_fasta = out_prefix + '_excluded.fasta'
        with open(good_fasta, 'w') as good, open(strange_fasta, 'w') as strange, open(singleton_fasta, 'w') as single,  open(excluded_fasta, 'w') as excluded:
            for rec in FastaIO.from_file(fasta):
                try:
                    cat, lineage = categorized[rec.id]
                    if cat == 'not_found':
                        print('Cluster ID not found in taxonomy: {} and not other taxonomy available'.format(rec.id))
                        continue

                    rec.id = '{};tax={};'.format(rec.id, utax_format(ranks, lineage, empty_patt))
                    if cat == 'good':
                        FastaIO.write(rec, good)
                    elif cat == 'singleton':
                        FastaIO.write(rec, single)
                    elif cat == 'excluded':
                        FastaIO.write(rec, excluded)
                    elif cat == 'strange':
                        FastaIO.write(rec, strange)
                    else:
                        print('error', file=sys.stderr)
                except KeyError:
                    print('not found: {}'.format(rec.id), file=sys.stderr)



p = argparse.ArgumentParser(description="""
This script tries to find reasonable taxonomic lineages for clusters of sequences with taxonomic annotation.
It takes a USEARCH/VSEARCH output file in UC format and a file containing the taxonomy, as produced
by extract_taxonomy.py. Output: '<prefix>_good.txt' with lineages that are probably correct.
strange: there are some contradictions within the taxonomy of the cluster. Singletons: only one sequence with taxonomic
annotations in cluster -> cannot know for sure if correct. excluded: Contains one of the '--exclude' keywords.
""")
p.add_argument("-o", "--out-prefix")
p.add_argument("-t", "--threshold", type=int, default=-2, help='Good lineages must be consistent below this rank number (negative numbers are counted from end)')
p.add_argument("-f", "--fasta", help="(optional) Fasta file with cluster centroids, FASTA files will be written for each group.")
p.add_argument("-e", "--exclude", default='environmental,uncultured')
p.add_argument("--empty-kws", default='Incertae')
p.add_argument("taxonomy", help="File with Genbank accession on each line")
p.add_argument("uc", type=argparse.FileType(), help="file in UCLUST format")
args = p.parse_args()

consensus_tax(**vars(args))
