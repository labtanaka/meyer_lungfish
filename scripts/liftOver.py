#!/usr/bin/env python3

import sys
import argparse
import gzip


class BEDFormatException(Exception):
    pass


def parseStructureFile(structfile):
    chrstruct = dict()
    hFile = None
    if structfile.endswith('.gz'):
        hFile = gzip.open(structfile,'rt')
    else:
        hFile = open(structfile)
    for strLine in hFile.readlines():
        scafID, contigID, strand, start, end = strLine.strip().split(' ')
        chrstruct[contigID] = {'scaffold': scafID, 'strand': strand, 'start': int(start), 'end': int(end) - 1}
    hFile.close()
    return (chrstruct)


"""
Converts the coordinates in the BED file. Since the orientation of the contigs within the scaffolds influences the
converted coordinates, throw an exception if the input BED file has less than 6 columns
"""
def convertBED(srcIn, chrstruct):
    with open(srcIn) as hFile:
        for strLine in hFile.readlines():
            arrCols = strLine.strip().split("\t")

            if len(arrCols) < 6:
                raise BEDFormatException("Cannot convert coordinates without the strand information")

            contigID = arrCols[0]
            iStart = int(arrCols[1])
            iEnd = int(arrCols[2])
            featStrand = arrCols[5]
            iThickStart = -1
            iThickEnd = -1
            if len(arrCols) >= 8:
                iStart = arrCols[6]
                iEnd = arrCols[7]

            # Find the scaffold that has the contig
            entry = chrstruct.get(contigID)
            if entry:
                # Convert the contig ID, the start and the end
                # If the strand of the corresponding contig is '+', simply convert the coordinates
                if entry['strand'] == '+':
                    iStart = entry['start'] + iStart - 1
                    iEnd = entry['start'] + iEnd - 1

                    if iThickStart > -1:
                        iThickStart = entry['start'] + iThickStart - 1
                        iThickEnd = entry['start'] + iThickEnd - 1
                # Otherwise, reverse the coordinates of the feature and change the strand
                else:
                    tmp = iStart
                    iStart = entry['end'] - iEnd + 1
                    iEnd = entry['end'] - tmp + 1
                    if featStrand == '+':
                        featStrand = '-'
                    else:
                        featStrand = '+'

                    if iThickStart > -1:
                        iThickStart = entry['start'] - iThickEnd + 1
                        iThickEnd = entry['start'] - iThickStart + 1

                contigID = entry['scaffold']

            else:
                print(f"WARN: the contig {contigID} is not contained within any scaffold. Left untouched", file=sys.stderr)
            
            # Print the new coordinates
            str = f"{contigID}\t{iStart}\t{iEnd}\t{arrCols[3]}\t{arrCols[4]}\t{featStrand}"
            if len(arrCols) >= 8:
                str += f"\t{iThickStart}\t{iThickEnd}"
            if len(arrCols) >= 9:
                str += "\t".join(arrCols[9,:])
            
            print(str, file=sys.stdout)

"""
Convert the coordinates in the GFF file.
"""
def convertGFF(srcIn, chrstruct):
    with open(srcIn) as hFile:
        for strLine in hFile.readlines():
            if strLine.startswith('#'):
                print(strLine, end='')
                continue
            
            contigID, progName, featClass, iStart, iEnd, score, featStrand, phase, attributes = strLine.strip().split("\t")
            # Find the scaffold that has the contig
            entry = chrstruct.get(contigID)
            if entry:
                if entry['strand'] == '+':
                    iStart = entry['start'] + int(iStart) - 1
                    iEnd = entry['start'] + int(iEnd) - 1
                # Otherwise, reverse the coordinates of the feature and change the strand
                else:
                    tmp = iStart
                    iStart = entry['end'] - int(iEnd) + 1
                    iEnd = entry['end'] - int(tmp) + 1
                    if featStrand == '+':
                        featStrand = '-'
                    else:
                        featStrand = '+'
                contigID = entry['scaffold']
            else:
                print(f"WARN: the contig {contigID} is not contained within any scaffold. Left untouched", file=sys.stderr)
            
            # Print the new coordinates
            print(f"{contigID}\t{progName}\t{featClass}\t{iStart}\t{iEnd}\t{score}\t{featStrand}\t{phase}\t{attributes}\t", file=sys.stdout)


"""
"""
def convertCoordinates(srcIn, chrstruct):
    with open(srcIn) as hFile:
        strLine = hFile.readline()
        if "gff-version" in strLine:
            print("The input file seems to be in GFF format", file=sys.stderr)
            convertGFF(srcIn, chrstruct)
        else:
            print("The input file seems to be in BED format", file=sys.stderr)
            convertBED(srcIn, chrstruct)


def main():
    ap = argparse.ArgumentParser(description='This script converts the coordinates in the contig assembly to those in the scaffolded assembly')
    ap.add_argument('--src', type=str, required=True, nargs=1, help='source coordinates in BED or in GFF3 format')
    ap.add_argument('--struct', type=str, required=True, nargs=1, help='path to the chromosome structure file')
    
    try:
        opts = ap.parse_args()
    except:
        ap.print_help()
        sys.exit(0)

    chrstruct = parseStructureFile(opts.struct[0])
    convertCoordinates(opts.src[0], chrstruct)

if __name__ == "__main__":
    main()