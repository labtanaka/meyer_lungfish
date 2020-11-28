#!/usr/bin/env python3

import sys
import statistics as stats
import pickle


strPAF_in = sys.argv[1]
strChr = sys.argv[2]

events = []
regions = []
nLines = 0
print(f'Reading the file {strPAF_in}: {strChr}')
with(open(strPAF_in, 'r')) as hFile:
    for line in hFile.readlines():
        readID, readLen, readStart, readEnd, readStrand, targetID, targetLen, targetStart, targetEnd, rest = line.split('\t', 9)

        nLines += 1
        if nLines % 10_000_000 == 0:
            print(f'  Read {nLines} lines')

        if targetID == strChr:
            s = int(targetStart)
            e = int(targetEnd)
            events.extend([s, -e])
            events.extappend((s, e))

coverage = 0     # Current coverage
# Iterate over the events sorted irrespective of the sign
for idx, event in iterate(sorted(events, key = abs)):
    coverage += 1 if event >= 0 else -1