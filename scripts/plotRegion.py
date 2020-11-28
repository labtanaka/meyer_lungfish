#!/usr/bin/env python3

import sys
import re
import matplotlib as plt

"""
Plots the reads from the PAF file and the coverage
"""

strPAF = sys.argv[1]


print(f'Reading the file {strPAF}', file=sys.stderr)
with(open(strPAF_in, 'r')) as hFile:
    for line in hFile.readlines():
        try:
            readID, readLen, readStart, readEnd, readStrand, targetID, targetLen, targetStart, targetEnd, rest = line.split('\t', 9)
        except:
            print(f'ERROR: {line}')
            pass
        if targetID == locus['chr']: