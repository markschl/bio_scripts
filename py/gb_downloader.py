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
p.add_argument("-i", "--ids", help="File with list of Genbank identifiers")
p.add_argument("-s", "--idsep", default=",", help="ID file separator (if multiple columns there)")
p.add_argument("-c", "--idcol", default=1, type=int, help="Column of file containing IDS")
p.add_argument("--blast", action="store_true", help="BLAST output ids")
p.add_argument("-q", "--query", help="Query ")
p.add_argument("-e", "--email", required=True, help="required...")
p.add_argument("-b", "--batch-size", type=int, default=20, help="batch size used to fetch data")
p.add_argument("-d", "--database", default="nucleotide", help="")
p.add_argument("-t", "--rettype", default="fasta", help="")
p.add_argument("-n", "--num-tries", default=20, type=int,
               help="Retry the given number of times before giving up if there is any error.")
p.add_argument("-o", "--outfile", type=argparse.FileType('w'), default='-')
args = p.parse_args()


Entrez.email = args.email


class Loader(object):
    def __init__(self, count, increment):
        self.count = count
        self.increment = increment
        self.position = 0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.next()

    def next(self):
        if self.position >= self.count:
            raise StopIteration()
        write_progress(
            f"downloading... {self.position:d} of {self.count:d} ({100 * self.position / self.count:.1f}%)")
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



if args.query:
    handle = Entrez.esearch(db=args.database, term=args.query, usehistory="y")
    props = Entrez.read(handle)
    handle.close()
    
    class SearchLoader(Loader):
        def fetch(self):
            return Entrez.efetch(
                db=args.database,
                webenv=props['WebEnv'], query_key=props["QueryKey"],
                rettype=args.rettype, retmode="text",
                retstart=self.position,
                retMax=args.batch_size
            )
            
    loader = SearchLoader(int(props['Count']), args.batch_size)

elif args.ids:
    f = open(args.ids)
    class IdLoader(Loader):
        def __init__(self, *arg):
            super(IdLoader, self).__init__(*arg)
            self.ids = None
            r = csv.reader(f, delimiter=args.idsep)
            idcol = args.idcol - 1
            self.iter = iter_chunks((self.get_id(row[idcol]) for row in r), self.increment)
            
        def get_id(self, id):
            return id
        
        def next(self):
            self.ids = list(next(self.iter))
            return super(IdLoader, self).next()
        
        def fetch(self):
            return Entrez.efetch(db=args.database, rettype=args.rettype, id=self.ids)
    loader = IdLoader

    if args.blast:
        class BlastIdLoader(IdLoader):
            def __init__(self, *arg):
                super(BlastIdLoader, self).__init__(*arg)
                self.pattern = re.compile(r"\w+\|(\d+)\|")

            def get_id(self, line):
                m = self.pattern.search(line)
                return m.group(1)


        loader = BlastIdLoader

    for i, line in enumerate(f): pass  # count lines
    f.seek(0)
    loader = loader(i + 1, args.batch_size)

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
            exit("\nDone\n")
        except Exception as  e:
            retry_in = 5 * (i + 1)
            print(f"\n{e}: retrying after {retry_in} seconds\n", file=sys.stderr, flush=True)
            time.sleep(retry_in)

    if i == args.num_tries - 1:
        print(f"Have been trying {args.num_tries} times, aborting and trying next. Fetching failed for items \
              {loader.position + 1} - {loader.position + args.batch_size + 1}", file=sys.stderr)
