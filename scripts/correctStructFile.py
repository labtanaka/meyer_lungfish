#!/usr/bin/env python3

import sys
import gzip


starts = {'NFORS_030035_pilon_pilon': 'scaf01.1',
          'NFORS_029962_pilon_pilon': 'scaf01.2',
          'NFORS_050262_pilon_pilon': 'scaf01.3',

          'NFORS_053044_pilon_pilon': 'scaf02.1',
          'NFORS_038576_pilon_pilon': 'scaf02.2',
          'NFORS_045284_pilon_pilon': 'scaf02.3',

          'NFORS_009231_pilon_pilon': 'scaf03.1',
          'NFORS_032622_pilon_pilon': 'scaf03.2',
          'NFORS_035245_pilon_pilon': 'scaf03.3',
          
          'NFORS_037603_pilon_pilon': 'scaf04.1',
          'NFORS_025729_pilon_pilon': 'scaf04.2'}

iOffset = 0
scafID = None
scafPrev = ''
with(gzip.open(sys.argv[1], 'rt')) as hFile:
    for line in hFile.readlines():
        fields = line.split()
        
        # If the current scaffold ID is different from the previous 
        # AND the current contig ID is present in the starts list,
        # then the scaffold must be split. Set the new scafID
        if starts.get(fields[1]):
            iOffset = int(fields[3]) - 1
            scafID = starts[fields[1]]
        elif fields[0] != scafPrev:
            iOffset = 0
            scafID = None
        scafPrev = fields[0]

        # If scafID not set, keep the scaffold ID as is and print the line.
        if not scafID:
            print(line, end='')
        else:
            # Otherwise, the scaffold ID must be changed and the coordinates adjusted.
            fields[0] = scafID
            fields[3] = str(int(fields[3]) - iOffset)
            fields[4] = str(int(fields[4]) - iOffset)
            print(" ".join(fields))

