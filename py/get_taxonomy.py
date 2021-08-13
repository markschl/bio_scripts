#!/usr/bin/env python3
import re
import sys
import csv
from collections import defaultdict
import xml.etree.ElementTree as ET
from urllib import request

from lib.util import iter_chunks


def get_taxonomy(source, output, db=None, accessions=None, file=None, file_delimiter='\t', missing_out=None, *args, **kwargs):
    if accessions is not None:
        accessions = accessions.split(',')
    elif file is not None:
        accessions = (row[0] for row in csv.reader(file, delimiter=file_delimiter))
    else:
        raise Exception('Either supply comma delimited accession list or file with accessions')

    writer = TaxWriter(output)

    if source == 'taxdb':
        missing = from_taxdb(db, accessions, writer, **kwargs)
    elif source == 'web':
        missing = from_web(accessions, writer, *args, **kwargs)
    else:
        raise Exception('Unknown source: {}'.format(source))

    if missing:
        if missing_out:
            missing_out.writelines(a + '\n' for a in missing)
        else:
            print('No taxonomy found for {} accessions. Use -m to obtain a list of them.'.format(len(missing)))


def from_taxdb(db, accessions, writer, **kw):
    from taxadb.accessionid import AccessionID

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
    from taxadb.schema import Accession, Taxa

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


def eutils(command, email, max_tries=20, **param):
    # ensure that there are no more than 3 requests per second
    t = time()
    t_diff = t - last_call['time']
    if t_diff < INTERVAL:
        sleep(INTERVAL - t_diff)
    last_call['time'] = t

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}.fcgi?email={}&tool={}'.format(command, email, TOOL)
    for k, v in param.items():
        url += '&{}={}'.format(k, v)
    #print(url, file=sys.stderr)
    for i in range(max_tries):
        try:
            response = request.urlopen(url)
        except Exception as e:
            retry_in = 5 * (i + 1)
            print(f"{e}: retrying after {retry_in} seconds\n", file=sys.stderr, flush=True)
            sleep(retry_in)
            if i == max_tries - 1:
                raise Exception(f"Too many attempts ({max_tries})")
    xml = response.read()
    return ET.fromstring(xml)


err_re = re.compile(r'Invalid uid ([^ ]+) at')


def from_web(accessions, writer, email=None, batch_size=100, **kw):
    if email is None:
        print('Email (-e) option is required with web source', file=sys.stderr)
        return

    no_taxid = []
    taxid2acc = defaultdict(list)
    for acc_chunk in iter_chunks(accessions, batch_size):
        acc_chunk = set(acc_chunk)
        result = eutils('esummary', email, db='nucleotide', id=','.join(acc_chunk), retmax=batch_size)
        for elem in result:
            if elem.tag == 'ERROR':
                acc = version = err_re.search(elem.text).group(1)
                assert acc is not None
                no_taxid.append(acc)
            else:
                acc = elem.find("Item[@Name='Caption']").text
                version = elem.find("Item[@Name='AccessionVersion']").text
                taxid = elem.find("Item[@Name='TaxId']").text
                # print('taxid', acc, elem.tag, elem.text, [(c.tag, c.text, c.get('Name')) for c in elem])
                taxid2acc[taxid].append(acc)
            try:
                acc_chunk.remove(acc)
            except KeyError:
                try:
                    acc_chunk.remove(version)
                except KeyError:
                    print(f'An unknown accession was returned from NCBI: {acc} (versioned: {version}).',
                          file=sys.stderr)

        if len(acc_chunk) == 0:
            print('Not all accessions returned, remaining: {}.'.format(','.join(acc_chunk), file=sys.stderr))

    for taxids in iter_chunks(taxid2acc.keys(), batch_size):
        taxids = set(taxids)
        result = eutils('efetch', email, db='taxonomy', id=','.join(taxids), retmax=batch_size)
        for elem in result.findall('Taxon'):
            taxid = elem.find('TaxId').text
            lineage = [
                (t.find('Rank').text, t.find('ScientificName').text)
                for t in elem.findall("./LineageEx/Taxon")
            ]
            lineage.append((elem.find('Rank').text, elem.find('ScientificName').text))
            for a in taxid2acc[taxid]:
                writer.write_lineage(a, lineage)
            taxids.remove(taxid)
        if len(taxids) == 0:
            print('Not all taxonomy IDs returned, remaining: {}.'.format(','.join(taxids), file=sys.stderr))
    return no_taxid


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser('''
    Script for retrieving taxonomic lineages for a set of Genbank accessions.
    Data is either retrieved from a Taxadb (https://github.com/HadrienG/taxadb)
    file in Sqlite format or using NCBI's Entrez webservice (for smaller sets).
    Using the Entrez webservice requires supplying a contact email address.
    ''')
    p.add_argument("-a", "--accessions", help="Comma delimited list of Genbank accessions")
    p.add_argument("-f", "--file",
                   help="File with Genbank accession on each line / in the first column if tab-delimited.",
                   type=argparse.FileType('r'), default='-')
    p.add_argument("--file-delimiter",
                   help='Delimiter for input file containing accessions ("--file"). Default: tab',
                   default='\t')
    p.add_argument("-o", "--output", type=argparse.FileType('w'), default='-')
    p.add_argument("-s", "--source", choices={'web', 'taxdb'}, default='web')
    p.add_argument("-d", "--db", help="Path to SQLITE database (taxdb source)")
    p.add_argument("-e", "--email", help="Required with web source")
    p.add_argument("-m", "--missing-out", type=argparse.FileType('w'))

    args = p.parse_args()

    get_taxonomy(**vars(args))
