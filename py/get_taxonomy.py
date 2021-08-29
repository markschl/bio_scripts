#!/usr/bin/env python3

import codecs
import gzip
import re
from sys import stderr
import csv
from collections import defaultdict
import xml.etree.ElementTree as ET
from urllib import request
from urllib.request import urlopen
from time import sleep, time

from remotezip import RemoteZip

from lib.util import iter_chunks


def get_taxonomy(source, output, db=None, accessions=None, file=None, file_delimiter='\t', missing_out=None, *args,
                 **kwargs):
    if accessions is not None:
        accessions = accessions.split(',')
    elif file is not None:
        accessions = [row[0] for row in csv.reader(file, delimiter=file_delimiter)]
    else:
        raise Exception('Either supply comma delimited accession list or file with accessions')

    writer = TaxWriter(output)

    if source == 'taxadb':
        missing = from_taxdb(db, accessions, writer, **kwargs)
    elif source == 'entrez':
        missing = from_entrez(accessions, writer, *args, **kwargs)
    elif source == 'ftp':
        missing = from_ftp(accessions, writer, *args, **kwargs)
    else:
        raise Exception('Unknown source: {}'.format(source))

    if missing:
        if missing_out:
            missing_out.writelines(a + '\n' for a in missing)
        else:
            print('No taxonomy found for {} accessions. Use -m to obtain a list of them.'.format(len(missing)))


def from_taxdb(db, accessions, writer, **kw):
    print("Fetching {} accessions from TaxaDB".format(len(accessions)), file=stderr)

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

last_call = {'time': 0}


def eutils(command, email, max_tries=100, **param):
    # ensure that there are no more than 3 requests per second
    t = time()
    t_diff = t - last_call['time']
    if t_diff < INTERVAL:
        sleep(INTERVAL - t_diff)
    last_call['time'] = t

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}.fcgi?email={}&tool={}'.format(command, email, TOOL)
    for k, v in param.items():
        url += '&{}={}'.format(k, v)
    # print(url, file=stderr)
    for i in range(max_tries):
        try:
            response = request.urlopen(url)
        except Exception as e:
            retry_in = 5 * (i + 1)
            print(f"{e}: retrying after {retry_in} seconds\n", file=stderr, flush=True)
            sleep(retry_in)
            if i == max_tries - 1:
                raise Exception(f"Too many attempts ({max_tries})")
            continue
        break
    xml = response.read()
    return ET.fromstring(xml)


err_re = re.compile(r'Invalid uid ([^ ]+) at')


def from_entrez(accessions, writer, email=None, id_batch_size=200, tax_batch_size=100, verbose=False, **kw):
    if email is None:
        print('Email (-e) option is required with web source', file=stderr)
        return

    if verbose:
        print("Fetching {} accessions from NCBI Entrez web service".format(len(accessions)), file=stderr)

    # First, obtain taxonomy IDs
    no_taxid = []
    taxid2acc = defaultdict(list)
    n = 0
    for acc_chunk in iter_chunks(accessions, id_batch_size):
        acc_chunk = set(acc_chunk)
        if verbose:
            n += len(acc_chunk)
            print('Obtaining taxonomy IDs ({} of {} = {:.1f} %)'.format(n, len(accessions), n / len(accessions) * 100),
                  file=stderr, end='\r')
        result = eutils('esummary', email, db='nucleotide', id=','.join(acc_chunk), retmax=id_batch_size)
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
                          file=stderr)

        if len(acc_chunk) > 0:
            print('Not all accessions returned in taxonomy ID search, remaining: {}.'.format(','.join(acc_chunk), file=stderr))

    # Obtain taxonomy from IDs
    n = 0
    print(taxid2acc)
    for taxids in iter_chunks(taxid2acc.keys(), tax_batch_size):
        taxids = set(taxids)
        if verbose:
            n += len(taxids)
            print('Obtaining lineages from IDs ({} of {} = {:.1f} %)'.format(n, len(accessions), n / len(accessions) * 100),
                  file=stderr, end='\r')
        result = eutils('efetch', email, db='taxonomy', id=','.join(taxids), retmax=tax_batch_size)
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
        if len(taxids) > 0:
            print('Not all taxonomy IDs returned in lineage search, remaining: {}.'.format(','.join(taxids), file=stderr))

    return no_taxid


def from_ftp(accessions, writer, accession_type='nucl_gb', verbose=False, **kw):
    if verbose:
        print("Fetching {} accessions from taxonomy dump on NCBI FTP server".format(len(accessions)), file=stderr)
    taxid2acc = defaultdict(list)
    accessions = {a: False for a in accessions}
    ftp_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy'
    bytes_stream = gzip.open(urlopen(f'{ftp_url}/accession2taxid/{accession_type}.accession2taxid.gz'))
    text_stream = codecs.iterdecode(bytes_stream, 'utf-8')
    n = 0
    for accession, version, taxid, _ in csv.reader(text_stream, delimiter='\t'):
        if accession in accessions:
            accessions[accession] = True
            taxid2acc[taxid].append(accession)
            n += 1
        if version in accessions:
            accessions[version] = True
            taxid2acc[taxid].append(version)
            n += 1
        if verbose and n % 100 == 0:
            print('Obtaining taxonomy IDs ({} of {} = {:.1f} %)'.format(n, len(accessions), n / len(accessions) * 100),
                  file=stderr, end='\r')

    # retrieve taxonomy
    # the ranks are not well documented
    ranks = ['no_rank', 'superkingdom', 'clade', 'kingdom', 'no_rank', 'phylum', 'subphylum',
             'class', 'subclass', 'order', 'suborder',
             'family', 'subfamily', 'genus', '?', '?', '?', 'species', 'varietas?']
    taxids = set(taxid2acc.keys())
    taxdump = RemoteZip(ftp_url + '/new_taxdump/new_taxdump.zip')
    byte_stream = taxdump.open('rankedlineage.dmp')
    text_stream = codecs.iterdecode(byte_stream, 'utf-8')
    n = 0
    for l in csv.reader(text_stream, delimiter='\t'):
        taxid = l[0]
        if taxid in taxids:
            n += 1
            if verbose and n % 100 == 0:
                print('Obtaining lineages ({} of {} = {:.1f} %)'.format(n, len(taxids),
                                                                            n / len(taxids) * 100),
                      file=stderr, end='\r')
            lineage = [(rank, name) for rank, name in zip(ranks, reversed(l)) if name != '' and name != '|']
            for acc in taxid2acc[taxid]:
                writer.write_lineage(acc, lineage)

    return [acc for acc, found in accessions.items() if not found]


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
    p.add_argument("-s", "--source", choices={'ftp', 'entrez', 'taxadb'}, default='ftp')
    p.add_argument("-t", "--accession-type",
                   choices={'nucl_gb', 'nucl_wgs', 'pdb', 'prot', 'dead_nucl', 'dead_prot', 'dead_wgs'},
                   default='nucl_gb')
    p.add_argument("-d", "--db", help="Path to SQLITE database (taxdb source)")
    p.add_argument("-e", "--email", help="Required with entrez source")
    p.add_argument("-m", "--missing-out", type=argparse.FileType('w'))
    p.add_argument("-v", "--verbose", action='store_true')

    args = p.parse_args()

    get_taxonomy(**vars(args))
