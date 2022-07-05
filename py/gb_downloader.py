#!/usr/bin/env python3

import argparse
import sys
import re
import csv
from lib.util import iter_chunks


from Bio import Entrez


p = argparse.ArgumentParser(
    description="This script allows downloading any data from Genbank, either by a search term (-q) or by specifying "
                "a file with Genbank identifiers or accessions (-i)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
p.add_argument("-a", "--accessions", type=argparse.FileType('r'), 
               help="File with list of Genbank accessions")
p.add_argument("-s", "--acc-sep", default=",", 
               help="accession file separator (if multiple columns present)")
p.add_argument("-c", "--acc-col", default=1, type=int, 
               help="Column of file containing accessions")
p.add_argument("-q", "--query", 
               help="Genbank Query to given database (nuccore by default, see -d)")
p.add_argument("-e", "--email", required=True, 
               help="E-mail address (required by Entrez service)")
p.add_argument("-b", "--batch-size", type=int, default=100, 
               help="batch size used to fetch data")
p.add_argument("-d", "--database", default="nuccore", 
               help="Database to query (default: nuccore)")
p.add_argument("-t", "--rettype", default="fasta", choices={"fasta", "gb", "acc", "seqid"},
               help="Return type (default: fasta); see https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and")
p.add_argument("-n", "--num-tries", default=100, type=int,
               help="Retry the given number of times before giving up if there is any error.")
p.add_argument("-o", "--outfile", type=argparse.FileType('w'), default='-')
args = p.parse_args()


Entrez.email = args.email


class Loader(object):
    def __init__(self, increment, count=None):
        self.count = count
        self.increment = increment
        self.position = 0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.next()

    def next(self):
        if self.count is None:
            write_progress(
                f"downloading... {self.position}")
        elif self.position <= self.count:
            write_progress(
                f"downloading... {self.position} of {self.count} ({100 * self.position / self.count:.1f}%)")
        data = self.fetch()
        self.position += self.increment
        return data

    def fetch(self):
        raise NotImplementedError()


def write_progress(progress, target=sys.stderr):
    """
    Allows writing progress messages to console.
    """
    target.write('\r')
    target.flush()
    target.write(str(progress))


class SearchLoader(Loader):
    def __init__(self, db, term, increment):
        handle = Entrez.esearch(db=db, term=term, usehistory="y")
        self.props = Entrez.read(handle)
        handle.close()
        super(SearchLoader, self).__init__(increment=increment, count=int(self.props['Count']))

    def fetch(self):
        if self.position >= self.count:
            raise StopIteration()
        return Entrez.efetch(
            db=args.database,
            webenv=self.props['WebEnv'], query_key=self.props["QueryKey"],
            rettype=args.rettype, retmode="text",
            retstart=self.position,
            retMax=args.batch_size
        )
 

class AccLoader(Loader):
    def __init__(self, accessions, *arg, **kwarg):
        super(AccLoader, self).__init__(*arg, **kwarg)
        self.ids = None
        self.iter = iter_chunks(accessions, self.increment)
        
    def get_id(self, id):
        return id
    
    def next(self):
        self.ids = list(next(self.iter))
        return super(AccLoader, self).next()
    
    def fetch(self):
        return Entrez.efetch(db=args.database, rettype=args.rettype, id=self.ids)


if args.query:    
    loader = SearchLoader(db=args.database, term=args.query, increment=args.batch_size)

elif args.accessions:
    f=args.accessions
    # first, count lines
    if f.seekable():
        for i, line in enumerate(f): pass
        f.seek(0)
        count = i + 1
    else:
        count = None
    # create accession iterator
    r = csv.reader(f, delimiter=args.acc_sep)
    c = args.acc_col - 1
    accessions = (row[c].strip() for row in r)
    accessions = (a for a in accessions if a)
    loader = AccLoader(accessions, increment=args.batch_size, count=count)

else:
    raise Exception("Please supply one of '-q' or '-a'")


import time

out = args.outfile

while True:
    # time.sleep(1)
    for i in range(args.num_tries):
        try:
            handle = next(loader) if i == 0 else loader.fetch()
            if handle is None:
                break
            data = handle.read()
            if not data.endswith("\n"):
                data += "\n"
            handle.close()
            out.write(data)
            break
        except StopIteration:
            print("\nDone\n", file=sys.stderr)
            exit()
        except Exception as  e:
            retry_in = 5 * (i + 1)
            print(f"\n{e}: retrying after {retry_in} seconds\n", file=sys.stderr, flush=True)
            time.sleep(retry_in)

    if i == args.num_tries - 1:
        print(f"Have been trying {args.num_tries} times, aborting and trying next. Fetching failed for items \
              {loader.position + 1} - {loader.position + args.batch_size + 1}", file=sys.stderr)
