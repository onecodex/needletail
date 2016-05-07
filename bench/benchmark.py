#!/usr/bin/env python
"""
Unscientific benchmarking of this versus the --release rust
implementation below using the %timeit Ipython magic (times in sec)

    n_kmers,  py_runtime, rust_runtime
    6594204,  14.4,       0.578

Both give identical counts on the files tested (and printing kmers out
and diff'ing the results gives no difference)
"""
from __future__ import print_function
import sys

from Bio.SeqIO import parse
from Bio.Seq import reverse_complement


def slid_win(seq, size=4, overlapping=True):
    """Returns a sliding window along self.seq."""
    itr = iter(seq)
    if overlapping:
        buf = ''
        for _ in range(size):
            buf += next(itr)
        for l in itr:
            yield buf
            buf = buf[1:] + l
        yield buf
    else:
        buf = ''
        for l in itr:
            buf += l
            if len(buf) == size:
                yield buf
                buf = ''

filename = sys.argv[1]

n_total = 0
n_canonical = 0

for s in parse(filename, 'fastq'):
    uppercase_seq = str(s.upper().seq)
    for kmer in slid_win(uppercase_seq, 4):
        canonical = min(kmer, reverse_complement(kmer))
        if canonical == 'CAGC':
            n_canonical += 1
        n_total += 1
        # print(canonical)

print(n_total, n_canonical)
