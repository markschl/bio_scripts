#!/usr/bin/env python3


import sys
import os
import re
import csv

import logging
import pysam
from lib.seq import reverse_complement, ambig
from lib import FastaIO
from lib.alignment import n_mismatches, compare, matches
from pprint import pprint

FNULL = open(os.devnull, 'w')

rev_ambig = {''.join(sorted(list(a))): b for b, a in ambig.items()}


def median(l):
    sorts = sorted(l)
    length = len(sorts)
    halflen = int(length / 2)
    if length % 2 == 0:
        return (sorts[halflen - 1] + sorts[halflen]) / 2
    return sorts[halflen]


def analyze(sam, reference, primers, output=sys.stdout, debug=None, skip=0, max_shift=3, f_key='fp', r_key='rp', **get_pos_args):
    if debug is not None:
        logger = logging.getLogger()
        logging.basicConfig(filename=debug, level=logging.DEBUG)
    else:
        logger = None

    fk = Key(f_key)
    rk = Key(r_key)

    w = csv.writer(output, delimiter='\t')
    # header
    w.writerow(('seqid', 'primer', 'direction', 'position', 'num_mismatches', 'mismatches',
                'num_read_mismatches', 'read_mismatches', 'mismatch_ratio', 'comment', 'is_problematic'))

    primer_map, region_info = parse_primerfile(primers)

    samfile = pysam.Samfile(sam)
    # get position of primers based on information from the sequence IDs
    # (which primer was used for sequencing)
    # This is done for all references found in the SAM file, regardless
    # of whether it is found in the reference FASTA
    seq_pdata = get_positions(samfile, primer_map, fk, rk, **get_pos_args)

    na = ''

    for rec in FastaIO.from_file(reference):
        seqid = rec.id

        if seqid not in seq_pdata:
            continue

        # if the FASTA matches the SAM references, there should be positional information
        positions = seq_pdata[seqid]

        n_problematic = len(positions)
        while True:
            # get mismatches of positioned primers
            for d in positions:
                d.set_mismatches(rec.seq, skip=0, max_shift=max_shift, logger=logger)
            # pprint(("mism", positions), sys.stderr)

            # predict additional primer positions based on existing positions
            # and known approx. positions of primers as defined in primer file
            # newly predicted pos. are returned
            pdata_predicted = predict_pos(positions, logger=logger)
            # pprint(("additional", pdata_predicted), sys.stderr)

            # get mismatches for the additional pos.
            for d in pdata_predicted:
                d.set_mismatches(rec.seq, skip=0, max_shift=max_shift, logger=logger)

            # check if placement of primers makes sense (domain order, known distance ranges)
            # misplaced primers are given a large position range
            problematic = validate(positions, region_info, logger=logger)
            # pprint(("problematic", problematic), sys.stderr)

            # get mismatches for those misplaced primers placed into a more confident range
            for d in problematic:
                d.set_mismatches(rec.seq, skip=0, max_shift=max_shift, logger=logger)

            problematic = [d for d in problematic if d.is_problematic()]
            # pprint(("end", positions), sys.stderr)

            # print "pr", len(problematic), n_problematic
            # if the number of problematic positions was reduced, then we need
            # to run the whole procedure again
            if len(problematic) == n_problematic:
                # no improvement to previous run
                break
            n_problematic = len(problematic)

        primer_positions = {p.primer.name: None if p.position is None else (p.position, p.position + p.primer.length) for p in positions}

        # write to output
        for d in positions:
            primer = d.primer
            d.problems = []

            if d.position is not None:
                set_read_mismatches(samfile, d, primer_positions, fk, rk)
                n_rm = sum(1 - (freqs['='] if '=' in freqs else 0) / sum(freqs.values())
                           for pos, freqs in d.read_mismatches)

                # ', '.join(d['problems'])
                w.writerow((seqid, primer.name, primer.get_direction(),
                            d.position + 1,
                            len(d.mismatches),
                            d.encode_mismatches(),
                            round(n_rm, 3),
                            d.encode_mismatch_freqs(),
                            d.mismatch_ratio,
                            'strange position' if d.strange_pos else na,
                            'x' if d.is_problematic() else na
                            ))
            elif d.checked:
                w.writerow((
                    seqid, primer.name, primer.get_direction(),
                    d.get_pos_or_guess(), na, na, na, na, na,
                    'strange position, approx. position' if d.strange_pos else 'approx. position',
                    'x' if d.is_problematic() else na
                ))
            else:
                w.writerow((seqid, primer.name, primer.get_direction(), na, na, na, na, na, na, na, na))


class Primer(object):
    # __slots__ = ('name', 'seq', 'is_reverse', 'length', 'region', 'offset')
    """
    Primer instances contain sequences...
    always on forward strand...

    region: name of conserved region
    offset: offset of primer in relation to region start (defined in header of input file)
    """

    def __init__(self, name, seq, region, offset, is_reverse=False):
        self.name = name
        self.seq = seq.split(',')
        self.is_reverse = is_reverse
        if is_reverse:
            self.seq = [reverse_complement(s) for s in self.seq]
        self.length = len(self.seq[0])
        self.region = region
        self.offset = offset

    def to_iupac(self):
        # mixes have linked bases, cannot be represented
        def get_iupac(variants):
            if len(variants) == 1:
                return variants[0]
            bases = ''.join(sorted(variants))
            if bases in rev_ambig:
                return rev_ambig[bases]
            return '?'

        return [get_iupac(b) for b in zip(*self.seq)]

    def get_direction(self):
        return 'r' if self.is_reverse else 'f'

    def align_seq(self, seq, logger=None):
        """
        Infers the exact position of primers by doing an ungapped alignment within
        the given range.
        """

        mers = [seq[i:i + self.length] for i in range(len(seq) - self.length + 1)]
        matches = [[compare(mer, pseq)[0] for pseq in self.seq] for mer in mers]
        n_mis_all = [[m.count(0) for m in match] for match in matches]
        best_i = [m.index(min(m)) for m in n_mis_all]
        n_mis = [n[i] for i, n in zip(best_i, n_mis_all)]

        best_n = min(n_mis)
        best_pos = n_mis.index(best_n)
        n_mis.remove(best_n)

        lh = len(mers) / 2
        # centering = abs(lh - best_pos) / lh
        mm_ratio = best_n / min(n_mis) if len(n_mis) > 0 and n_mis[0] != 0 else 0

        # if best_n >= 6 and logger is not None:
        #     # debug
        #     logger.debug("\n{}: mismatch ratio {:.2f}, {} {}\n{} seq\n{} primer".format(
        #         self.name,
        #         abs(lh - best_pos) / lh,
        #         mm_ratio,
        #         (abs(lh - best_pos) / lh) * mm_ratio,
        #         seq,
        #         ','.join(self.seq),
        #         '\n'.join('{}: {}'.format(mer, min(n_mismatches(mer, s) for s in self.seq)) for mer in mers)
        #     ))

        best_mer = seq[best_pos:best_pos + self.length]
        best_matches = matches[best_pos][best_i[best_pos]]
        mismatches = [(i, best_mer[i]) for i, m in enumerate(best_matches) if not m]

        return best_n, best_pos, mismatches, mm_ratio


class PrimerPosition(object):
    # __slots__ = ('ref', 'primer', 'problems', 'checked', 'read_mismatches', 'read_range', 'mismatches', 'num_mismatches', 'mismatch_ratio', 'strange_pos', '_pos_guess', 'position')
    """
    Primer position on a specific sequence
    """

    def __init__(self, reference_name, primer):
        self.ref = reference_name
        self.primer = primer
        self.problems = []
        self.checked = False
        self.read_mismatches = None
        self.read_range = None
        self.mismatches = None
        self.num_mismatches = None
        self.mismatch_ratio = None
        self.strange_pos = None
        self._pos_guess = None
        self.position = None

    def position_known(self, start, end):
        self.checked = True
        self._pos_guess = (start, end)

    def get_pos_or_guess(self):
        if self.position is None:
            return self._pos_guess[round((len(self._pos_guess) - 1) / 2)]
        return self.position

    def set_mismatches(self, seq, max_shift=3, skip=0, logger=None):
        """
        Sets mismatches, and based on the number of mismatches at 'max_shift'
        adjacent positions, the guessed positions are confirmed.
        Sets the num_mismatches, mismatches, mismatch_ratio and position attributes

        Args:
            seq:
            max_shift:
            skip:

        Returns: dict { primer_name: (num_mismatches, mismatch_code, 1-based_position), ... }

        """

        if not self.checked:
            return False

        start, end = self._pos_guess

        # start and end of slice
        start = max(start - max_shift, skip - 1)  # limit at start
        end = min(len(seq) - skip, end + max_shift + self.primer.length)  # limit at end
        subseq = seq[start:end]
        # TODO: the subseq should actually never be empty if end-start > 0, not sure why this happens

        if len(subseq) < self.primer.length:  # pos0 < skip or pos0 + plen + 1 > seqlen - skip or end - start <= 0:
            #self.checked = False
            return False

        self.num_mismatches, offset, self.mismatches, self.mismatch_ratio = self.primer.align_seq(subseq, logger=logger)
        self.position = start + offset

    def is_problematic(self):
        if self.num_mismatches is None:
            # no information -> assume the best
            return False
        return self.num_mismatches > 3 and self.mismatch_ratio > 0.5 or self.strange_pos == True

    def encode_mismatches(self):
        primer_seq = self.primer.to_iupac()

        # assemble: 1. position, counted from end; 2: primer base; 3: sorted base frequencies
        if not self.primer.is_reverse:
            m = [(pos, primer_seq[pos], base) for pos, base in self.mismatches]
        else:
            m = reversed([(self.primer.length - pos - 1, reverse_complement(primer_seq[pos]), reverse_complement(base))
                          for pos, base in self.mismatches])

        return ';'.join('{}:{}>{}'.format(pos + 1, ref_base, base)
                        for pos, ref_base, base in m)

    def encode_mismatch_freqs(self, pos_delim=':', ref_delim='>'):
        primer_seq = self.primer.to_iupac()

        # assemble: 1. position, counted from end; 2: primer base; 3: sorted base frequencies
        if not self.primer.is_reverse:
            m = [(pos, primer_seq[pos],
                  sorted(list(base_freqs.items()), key=lambda x: x[1], reverse=True))
                 for pos, base_freqs in self.read_mismatches]
        else:
            # reverse positions again, and reverse complement sequences
            m = [(self.primer.length - pos - 1, reverse_complement(primer_seq[pos]),
                 [(reverse_complement(b), f) for b, f in sorted(list(base_freqs.items()), key=lambda x: x[1], reverse=True)])
                 for pos, base_freqs in self.read_mismatches]

        return ';'.join('{}:{}>{}'.format(
            pos + 1,
            ref_base,
            ','.join('{}{}'.format(count, base) for base, count in base_freqs)
        ) for pos, ref_base, base_freqs in m)


class Key(object):
    def __init__(self, key):
        self.patt = re.compile("%s=([^;$]+)" % re.escape(key))

    def get_key(self, string):
        m = self.patt.search(string)
        if m is None:
            return ''
        return m.group(1)


def validate(positions, region_info, logger=None):
    """
    Check if placement of primers makes sense (domain order, known distance ranges)
    misplaced primers are given a large position range.

    """
    # { region_name: [PrimerPosition, ...] }
    by_region = {}
    for d in positions:
        if d.checked:
            reg = d.primer.region
            if reg not in by_region:
                by_region[reg] = []
            by_region[reg].append(d)

    if not by_region:
        return []

    region_list = [(rname, rinfo, by_region[rname]) for rname, rinfo in region_info if rname in by_region]
    scores = [iter_mean(d.mismatch_ratio if d.mismatch_ratio is not None else 0 for d in pdat) - len(pdat) * 0.01 for rname, rinfo, pdat in region_list]
    anchor = scores.index(min(scores))

    rname0, rng0, pdat0 = region_list[anchor]
    # approx. left bound of region (mean from primer systems) in sequence
    left_bound0 = int(iter_mean(d.get_pos_or_guess() - d.primer.offset for d in pdat0))

    problematic = []

    # go back sequentially
    offset = [0, 0]  # cumulative left offset of region compared to anchor region
    for rname, rng, pdat in reversed(region_list[:anchor]):
        diffs = []
        offset[0] += rng[0]
        offset[1] += rng[1]

        for d in pdat:
            off = d.primer.offset
            left_bound = d.get_pos_or_guess() - off
            diff = left_bound0 - left_bound

            if offset[0] <= diff <= offset[1]:
                diffs.append(diff)
                d.strange_pos = False
            else:
                if logger is not None:
                    logger.debug("%s: strange pos. of %s in %s: %d should be in %d-%d" % (
                        d.ref, d.primer.name, rname, diff, offset[0], offset[1]))

                o = left_bound + off
                d.position_known(o - offset[0], max(0, o - offset[1]))
                d.strange_pos = True
                problematic.append(d)
        if diffs:
            # there are some primers that seem to be ok, therefore we can set the offset
            # with some confidence
            diff = int(iter_mean(diffs))
            offset = [diff, diff]

    if rng0 is not None:
        # go forward sequentially
        offset = list(rng0)  # cumulative offset range
        for rname, rng, pdat in region_list[anchor + 1:]:
            diffs = []
            for d in pdat:
                off = d.primer.offset
                left_bound = d.get_pos_or_guess() - off
                diff = left_bound - left_bound0
                if offset[0] <= diff <= offset[1]:
                    diffs.append(diff)
                    d.strange_pos = False
                else:
                    if logger is not None:
                        logger.debug("%s: strange pos. of %s in %s: %d should be in %d-%d" % (
                            d.ref, d.primer.name, rname, diff, offset[0], offset[1]))
                    o = left_bound0 + off
                    d.pos_guess = [offset[0] + o, offset[1] + o]
                    d.strange_pos = True
                    problematic.append(d)
            if rng is not None:
                # regions following, add current region to offset range
                if diffs:
                    # there are some primers that seem to be ok, therefore we can set the offset
                    # with some confidence
                    diff = int(iter_mean(diffs))
                    offset[0] += diff
                    offset[1] += diff
                else:
                    offset[0] += rng[0]
                    offset[1] += rng[1]

    return problematic


def iter_mean(x):
    """
    Mean value from any iterable (not just lists)

    Args:
        x: iterable yielding numeric values

    Returns: float

    """
    n = 0
    s = 0
    for v in x:
        s += v
        n += 1
    return s / n


def parse_primerfile(primerfile, logger=None):
    with open(primerfile) as p:
        line1 = next(p).strip()
        assert line1.startswith('#'), 'Expecting region information at start, preceded by #'
        primers = list(csv_dictreader(p, delimiter='\t'))

    # parse info on regions
    info = iter(line1[1:].lstrip().split('\t'))
    region_info = []
    while True:
        try:
            region_name = next(info)
            spacer_info = next(info).split('+')
            region_end_rng = [int(x.strip()) for x in spacer_info[0].split('-', 1)]
            spacer_rng = [int(x.strip()) for x in spacer_info[1].split('-', 1)]
            region_info.append((region_name, (
                min(region_end_rng) + min(spacer_rng),
                max(region_end_rng) + max(spacer_rng)
            )))

        except StopIteration:
            region_info.append((region_name, None))
            break

    if logger is not None:
        logger.debug('region info: {} -> {}'.format(line1, region_info))

    # Primer instances

    primer_map = {
        p['name']: Primer(
            p['name'],
            p['seq'],
            p['region'],
            int(p['offset']),
            is_reverse=p['direction'] == 'r'
        ) for p in primers
        }

    return primer_map, region_info


def csv_dictreader(f, **kwargs):
    header = next(csv.reader(f, **kwargs))
    for r in csv.DictReader(f, header, **kwargs):
        yield {k: r[k] for k in header}


def get_positions(samfile, primer_map, f_key, r_key, freq_threshold=0.2):
    """
    Parses forward and reverse primer names of each SAM reference sequence
    -> list of initialized PrimerPosition objects
    Args:
        samfile:
        primer_map: dict { 'primer name': Primer, ... }
        f_key: Key instance
        r_key: Key instance
        freq_threshold: float

    Returns: dict { reference_name: [ PrimerPosition ] }

    """
    # get start and stop positions of the forward and reverse primers
    # that are stored in the sequence ID (retrieved using Key.get_key())
    # { reference_name: { primer_name: [ { pos: num_occurrences }, direction ('f'/'r') ] } }
    data = {}
    for read in samfile:
        if not read.is_unmapped:
            bl = read.get_blocks()
            start = bl[0][0]
            stop = bl[-1][1]

            c = read.cigartuples
            if c[0][0] == 1:
                # insertion at start -> subtract
                # TODO: correct? works with VSEARCH
                start -= c[0][1]

            fprimer = f_key.get_key(read.query_name)
            rprimer = r_key.get_key(read.query_name)

            ref = read.reference_name
            if ref not in data:
                data[ref] = {}
            positions = data[ref]

            for primer_name, pos, is_reverse in ((fprimer, start, False), (rprimer, stop, True)):
                if primer_name not in positions:
                    positions[primer_name] = [{}, is_reverse]
                primer_pos, primer_is_rev = positions[primer_name]

                assert primer_is_rev == is_reverse, "Error: Primer with name '{}' was found as forward AND reverse primer".format(
                    primer_name)

                if pos not in primer_pos:
                    primer_pos[pos] = 0
                primer_pos[pos] += 1

    out = {}
    for ref, positions in data.items():
        # positions: dict { primer_name: [ { pos: num_occurrences }, direction ] }
        p_out = out[ref] = []
        for primer_name, primer in primer_map.items():
            d = PrimerPosition(ref, primer)
            p_out.append(d)
            if primer_name in positions:
                pos, is_reverse = positions[primer_name]
                assert is_reverse == primer.is_reverse
                (start, end) = get_pos_rng(pos, freq_threshold)
                if not is_reverse:
                    # positions are always in relation to forward strand
                    start -= primer.length
                    end -= primer.length
                d.position_known(start, end)
    return out


def get_pos_rng(pos, freq_threshold=0.2):
    """
    Returns a range of possible start/end positions of primers based on a dict
    containing actual found positions
    Args:
        pos: dict { position: frequency }
        freq_threshold: min. frequency required in order to be in range

    Returns: start, end (1-based?)

    """
    pos = sorted(list(pos.items()), key=lambda x: x[1], reverse=True)
    best_pos, best_count = pos[0]
    outpos = {best_pos}
    for p, count in pos[1:]:
        if count / best_count >= freq_threshold:
            outpos.add(p)
    outpos = sorted(list(outpos))
    return outpos[0], outpos[-1]


def primers_by_region(primers):
    out = {}
    for p in primers:
        region = p.region
        if region not in out:
            out[region] = []
        out[region].append(p)
    return out


def predict_pos(positions, max_diff=5, logger=None):
    by_region = {}
    for d in positions:
        region = d.primer.region
        if region not in by_region:
            by_region[region] = []
        by_region[region].append(d)

    predicted = []

    for region, pdat in by_region.items():
        sorted_mis = []
        # make list of primer systems sorted by number of mismatches.
        # those not in the mismatches dict are added as 'new' with their positions
        # guessed from the offset value
        # positions with many mismatches that do not fit with the offset are returned as 'new' as well
        not_found = []

        for d in pdat:
            if d.checked:
                sorted_mis.append((d.num_mismatches, d))
            else:
                not_found.append(d)

        sorted_mis.sort(key=lambda m: m[0] if m[0] is not None else 0)

        if sorted_mis:
            d0 = sorted_mis[0][1]
            mean_diff = d0.get_pos_or_guess() - d0.primer.offset
            for _, d in sorted_mis:
                diff = d.get_pos_or_guess() - d.primer.offset
                if abs(diff - mean_diff) > max_diff:
                    if logger is not None:
                        logger.debug('strange position: {}, {}, {}'.format(d.ref, d.position, d.primer.offset))
                    not_found.append(d)
                else:
                    mean_diff = (mean_diff + diff) // 2

            for d in not_found:
                pos = mean_diff + d.primer.offset
                d.position_known(pos, pos)
                predicted.append(d)

    return predicted


def set_read_mismatches(samfile, d, primer_coords, f_key, r_key, freq_threshold=0.01):
    """

    """
    if d.position is None:
        # set_mismatches() found that there is no sequence information
        # for this primer
        return

    positions = []
    mismatches = []

    primer = d.primer
    pseqs = list(zip(*primer.seq))
    start = d.position
    end = start + primer.length

    def in_for_primer_range(pos, id):
        ppos = primer_coords[f_key.get_key(id)]
        # print('for', id, f_key.get_key(id), pos, ppos, end='')
        if ppos is not None:
            return pos < ppos[1]
        return False

    def in_rev_primer_range(pos, id):
        ppos = primer_coords[r_key.get_key(id)]
        # print('rev', id, r_key.get_key(id), pos, ppos, end='')
        if ppos is not None:
            return pos >= ppos[0]
        return False

    in_primer_range = in_rev_primer_range if d.primer.is_reverse else in_for_primer_range

    for column in samfile.pileup(d.ref, start, end):
        pos = column.reference_pos
        offset = pos - start

        if start <= pos < end:
            positions.append(pos)
            base_map = {}
            n = 0
            for p in column.pileups:
                if p.query_position is None or in_primer_range(p.alignment.reference_start + p.query_position, p.alignment.query_name):
                    # seq may be unaligned (-> p.query_position = None)
                    # or there may be some 'incorrect' bases in primer range, that should be removed
                    continue
                n += 1
                qpos = p.query_position
                if qpos is not None:
                    b = p.alignment.query_sequence[qpos]
                    if b in base_map:
                        base_map[b] += 1
                    else:
                        base_map[b] = 1

            mm_count = {}
            n_mismatches = 0
            for base, count in base_map.items():
                if any(matches(base, b) for b in pseqs[offset]) or count / n < freq_threshold:
                    if '=' not in mm_count:
                        mm_count['='] = 0
                    mm_count['='] += count
                else:
                    mm_count[base] = count
                    n_mismatches += 1

            if n_mismatches > 0:
                mismatches.append((offset, mm_count))

    d.read_mismatches = mismatches
    d.read_range = (positions[0], positions[-1]) if positions else None


if __name__ == '__main__':
    import argparse

    p = argparse.ArgumentParser(description="""
    Infers the number and positions of primer mismatches from a BAM/SAM alignment.
    The aligned sequence (query) IDs should contain the primer names in a
    key=name format (by default, fp is assumed for forward primer and rp for rev.
    primer, change with -f/-r). In addition, the primers should be removed, the
    sequence should start/end at correct positions. If the ends are fuzzy,
    increase the search range using --max-shift. The best positions are inferred
    based on the mismatch distribution (there should be an optimal position).

    Additional certainty is given by the primer file, which specifies how the
    positions of the primers relative to a conserved region should be. The header
    line (starting with a #) contains the order and position of these regions.
    The format is: <region_name> <region_length> + <min_dist> - <max_dist> <region2_name>...
    (min_dis/max_dist is a range of possible inter-region distances). If regions
    are in the wrong order or the distances don't fit, there will be a 'x' in the
    strange_pos column of the output.

    Example primer file:
    # SSU	23 + 70-420	5.8S	157 + 100-400	LSU
    name	seq	direction	region	offset
    58S_left	RRCAAYGGATCTCTHGG	f	5.8S	9
    (...)

    Output columns:
    seqid: reference ID (should match in reference FASTA and SAM/BAM file)
    primer: primer name
    direction: f=forward, r=reverse
    position: inferred 1-based position relative to SAM reference (FASTA reference?)
    num_mismatches: number of mismatches to reference
    mismatches: Encoded mismatch positions
    num_read_mismatches: number of mismatches inferred from all aligned sequences
    read_mismatches: Encoded distribution of mismatches in aligned sequences
    mismatch_ratio: ratio best_mismatch / median(mismatches)...
    strange_pos: 'x' if the inferred position seems not to fit the domain configuration / distances
    is_problematic: 'x' if many mismatches and therefore not sure, if the placement is correct
    """)
    p.add_argument("sam", help="SAM or BAM file")
    p.add_argument("reference", help="FASTA file with ref. sequences.")
    p.add_argument("primers", help="Primer file.")
    p.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('wb'), help="defaults to stdout")
    p.add_argument("-d", "--debug",
                   help="Debug output file (default: none)")
    p.add_argument("-t", "--freq_threshold", default=0.2, type=float,
                   help="Include sequence positions with at least the given frequency in the output.")
    p.add_argument("-f", "--f_key", default='fp', help="Key that specifies the forward primer name (default: 'fp')")
    p.add_argument("-r", "--r_key", default='rp', help="Key that specifies the reverse primer name (default: 'rp')")
    p.add_argument("-s", "--max_shift", type=int, default=3,
                   help="Number of positions to test left and right of primer (default: 3)")
    p.add_argument("-k", "--skip", type=int, default=0,
                   help="Number of bases to skip at ends of reference sequences (may be primers)")
    args = p.parse_args()

    analyze(**vars(args))
