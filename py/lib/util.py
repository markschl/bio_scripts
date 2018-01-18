
from itertools import islice, chain


def iter_chunks(iterable, size):
    """Splits an iterable into even sized chunks"""
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))


def utax_format(ranks, lineage, empty_patt):
    return ','.join('{}:{}'.format(r[0], name.replace(',', '_'))
                    for r, name in zip(ranks, lineage)
                    if len(clean_name(name, empty_patt)) > 0)


def clean_name(name, patt):
    if patt.search(name) is not None:
        return ''
    return name
