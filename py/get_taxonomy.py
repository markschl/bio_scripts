#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict
import xml.etree.ElementTree as ET
from urllib import request

from taxadb.schema import Accession, Taxa
from taxadb.accessionid import AccessionID

from lib.util import iter_chunks


def tax(source, output, db=None, accessions=None, file=None, missing_out=None, *args, **kwargs):

    if accessions is not None:
        accessions = accessions.split(',')
    elif file is not None:
        accessions = (line.rstrip('\r\n ') for line in file if len(line) >= 2)
    else:
        raise Exception('Either supply comma delimited accession list or file with accessions')

    writer = TaxWriter(output)

    if source == 'taxdb':
        missing = from_taxdb(db, accessions, writer, **kwargs)
    elif source == 'web':
        missing = from_web(accessions, writer, *args, **kwargs)
    else:
        print('unknown source:', source)

    if missing_out and missing:
        missing_out.writelines(a + '\n' for a in missing)


def from_taxdb(db, accessions, writer, **kw):

    accdb = AccessionID(dbtype='sqlite', dbname=db)

    accessions = (a.split('.')[0] for a in accessions)
    missing = []
    
    for acc in iter_chunks(accessions, 500):
        acc = list(acc)
        acc_set = set(acc)
        for a, lineage in lineage_names_levels(accdb, acc):
            acc_set.remove(a)
            writer.write_lineage(a, reversed(lineage))
        missing += list(acc_set)
    return missing
    

def lineage_names_levels(accdb, acc_number_list):
    """Get a lineage name for accession ids

    Given a list of acession numbers, yield the accession number and their
        associated lineage as tuples

    Args:
        acc_number_list (:obj:`list`): a list of accession numbers

    Yields:
        tuple: (accession id, lineage name)

    """
    accdb.check_list_ids(acc_number_list)
    with accdb.db.atomic():
        query = Accession.select().where(
            Accession.accession << acc_number_list)
        for i in query:
            lineage_list = []
            current_lineage = (i.taxid.lineage_level, i.taxid.tax_name)
            parent = i.taxid.parent_taxid
            while current_lineage[1] != 'root':
                lineage_list.append(current_lineage)
                new_query = Taxa.get(Taxa.ncbi_taxid == parent)
                current_lineage = (new_query.lineage_level, new_query.tax_name)
                parent = new_query.parent_taxid
            yield (i.accession, lineage_list)


class TaxWriter(object):
    def __init__(self, f):
        self.writer = csv.writer(f, delimiter='\t')
    
    def write_lineage(self, accession, lineage):
        self.writer.writerow([accession] + ['{}:{}'.format(rank.replace(' ', '_'), name) for rank, name in lineage])    


TOOL = 'get_taxonomy.py'
INTERVAL = 0.35

from time import sleep, time

last_call = {'time': 0}
def eutils(command, email, **param):
    # ensure that there are no more than 3 requests per second
    t = time()
    t_diff = t - last_call['time']
    if t_diff < INTERVAL:
        sleep(INTERVAL - t_diff)
    last_call['time'] = t
    
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}.fcgi?email={}&tool={}'.format(command, email, TOOL)
    for k, v in param.items():
        url += '&{}={}'.format(k, v)
    response = request.urlopen(url)
    xml = response.read()
    return ET.fromstring(xml)


def from_web(accessions, writer, email=None, **kw):

    if email is None:
        print('Email (-e) option is required with web source', file=sys.stderr)
        return

    no_taxid = []
    taxid2acc = defaultdict(list)
    
    for acc in iter_chunks(accessions, 200):
        acc = list(set(acc))
        result = eutils('esummary', email, db='nucleotide', id=','.join(acc), retmax=200)

        for a, elem in zip(acc, result.getchildren()):
            if elem.tag == 'ERROR':
                no_taxid.append(a)
            else:
                taxid = elem.find("Item[@Name='TaxId']").text
                taxid2acc[taxid].append(a)

    for taxids in iter_chunks(taxid2acc.keys(), 200):
        taxids = list(taxids)
        result = eutils('efetch', email, db='taxonomy', id=','.join(taxids), retmax=200)
        for taxid, tax in zip(taxids, result.findall('Taxon')):
            lineage = [
                (t.find('Rank').text, t.find('ScientificName').text)
                for t in tax.findall("./LineageEx/Taxon")
            ]
            lineage.append((tax.find('Rank').text, tax.find('ScientificName').text))

            for a in taxid2acc[taxid]:
               writer.write_lineage(a, lineage)

    return no_taxid



import argparse


p = argparse.ArgumentParser('''
Script for retrieving taxonomic lineages for a set of Genbank accessions.
Data is either retrieved from a Taxadb (https://github.com/HadrienG/taxadb)
file in Sqlite format or using NCBI's Entrez webservice (for smaller sets)
''')
p.add_argument("-a", "--accessions", help="Comma delimited list of Genbank accessions")
p.add_argument("-f", "--file", help="File with Genbank accession on each line",
               type=argparse.FileType('r'), default='-')
p.add_argument("-o", "--output", type=argparse.FileType('w'), default='-')
p.add_argument("-s", "--source", choices={'taxdb', 'web'}, default='taxdb')
p.add_argument("-d", "--db", help="Path to SQLITE database (taxdb source)")
p.add_argument("-e", "--email", help="Required with web source")
p.add_argument("-m", "--missing-out", type=argparse.FileType('w'))

args = p.parse_args()

tax(**vars(args))
