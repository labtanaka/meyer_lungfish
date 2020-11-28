#!/usr/bin/env python3

import sys
import statistics as stats
import pickle
import array

strPAF_in = sys.argv[1]
strChr = sys.argv[2]

cov_data = None
nLines = 0
print(f'Reading the file {strPAF_in}: {strChr}')
with(open(strPAF_in, 'r')) as hFile:
    for line in hFile.readlines():
        readID, readLen, readStart, readEnd, readStrand, targetID, targetLen, targetStart, targetEnd, rest = line.split('\t', 9)

        nLines += 1
        if nLines % 10_000_000 == 0:
            print(f'  Read {nLines} lines')

        if targetID == strChr:
            if not cov_data:
                cov_data = [0] * int(targetLen)
            for i in range(int(targetStart), int(targetEnd)):
                cov_data[i] += 1

print(f'Calculating statistics')
mn = stats.mean(cov_data)
sd = stats.stdev(cov_data)
print(f'Mean coverage: {mn}')
print(f'STD dev:       {sd}')

with(open(f'{strChr}.coverage.pkl', 'wb')) as hOut:
    pickle.dump(cov_data, hOut)