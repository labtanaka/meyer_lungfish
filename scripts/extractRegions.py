#!/usr/bin/env python3

import sys
import re

strPAF_in = sys.argv[1]
strPAF_out = sys.argv[2]
strLocus = sys.argv[3]
strPDF_out = sys.argv[4]

locus = {'chr': strLocus, 'start': 0, 'end': -1}

m = re.search('([^:]+):([0-9]+)-([0-9]+)', strLocus)
if m:
    locus['chr'] = m.group(1)
    locus['start'] = int(m.group(2))
    locus['end'] = int(m.group(3))
print(f"Extracting the data for {locus['chr']}:{locus['start']}-{locus['end']}", file=sys.stderr)

print(f'Reading the file {strPAF_in}', file=sys.stderr)
lines = []
coord_lims = [-1, -1]
with(open(strPAF_in, 'r')) as hFile:
    for line in hFile.readlines():
        try:
            readID, readLen, readStart, readEnd, readStrand, targetID, targetLen, targetStart, targetEnd, rest = line.split('\t', 9)
        except:
            print(f'ERROR: {line}')
            pass
        if targetID == locus['chr']:
            if locus['end'] == -1:
                locus['end'] = int(targetLen)

            if coord_lims[0] == -1:
                coord_lims[0] = int(targetStart)
                coord_lims[1] = int(targetEnd)
             
            if locus['start'] <= int(targetStart) <= locus['end'] or locus['start'] <= int(targetEnd) <= locus['end']:
                lines.append([readID, readLen, readStart, readEnd, readStrand, targetID, targetLen, targetStart, targetEnd, rest])
                if int(targetStart) < coord_lims[0]:
                    coord_lims[0] = int(targetStart)
                if int(targetEnd) > coord_lims[1]:
                    coord_lims[1] = int(targetEnd)

print(coord_lims, file=sys.stderr)
print(f'Adjusting the coordinates', file=sys.stderr)
for line in lines:
    line[6] = str(coord_lims[1] - coord_lims[0])
    line[7] = str(int(line[7]) - coord_lims[0])
    line[8] = str(int(line[8]) - coord_lims[0])
    print("\t".join(line), end='')

