#!/usr/bin/env python3

import matplotlib.pyplot as plt
import statistics as stats
import functools
import sys
import pickle
import math
from matplotlib.lines import Line2D
import numpy as np
from matplotlib import gridspec
import pandas as pd

strPkl = sys.argv[1]
strScaffoldID = sys.argv[2]
nBinSize = int(sys.argv[3])
strRepeatMasker = sys.argv[4]  #'/groups/tanaka/Projects/lungfish-genome/manuscript/Revision_II/NFOR_v3.1.repeatmasker.SCAFFOLDS.gff3'
strOutputPDF = sys.argv[5]

print("PARAMETERS:")
print(f'Pickle file: {strPkl}')
print(f'ScaffoldID: "{strScaffoldID}"')
print(f'Bin size: {nBinSize}')
print(f'RepeatMasker file: {strRepeatMasker}', flush=True)

##############
def compare(r1, r2):
    if r1['start'] < r2['start']:
        return -1
    if r1['start'] > r2['start']:
        return 1
    if r1['end'] < r2['end']:
        return -1
    if r1['end'] > r2['end']:
        return 1
    return 0
###############


print(f'Loading the file {strPkl}')
dat = None
with(open(strPkl, 'rb')) as hFile:
    dat = pickle.load(hFile)

print(f'Loaded {len(dat):,} items')

print('Calculating the mean and the standard deviation')
mn = stats.mean(dat)
sd = stats.stdev(dat)
print(f'Mean: {mn}')
print(f'SD: {sd}')


print('Binning the data')
bins = [0] * math.ceil(len(dat) / nBinSize) 
iBin = 0
iPos = 0
while iPos < len(dat):
    bindata = dat[iPos : iPos + nBinSize]
    iPos += nBinSize
    bins[iBin] = stats.mean(bindata)
    iBin += 1
    
    if iBin % 1000 == 0:
        print(f'  Processed {iBin} bins')

print(f'Plotting the data')
bin_mean = stats.mean(bins)
bin_std = stats.stdev(bins)

fig = plt.figure(figsize=(20,15))
plt.rc('font', family='sans-serif') 
plt.rc('font', serif='Helvetica') 
plt.rc('text', usetex='false')
#gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1])
#axs = [plt.subplot(gs[0]), plt.subplot(gs[1])] 
gs = gridspec.GridSpec(1, 1)
axs = [plt.subplot(gs[0])]

bins_log = [math.log(x) if x > 0 else 0 for x in bins]
axs[0].fill_between(range(0, len(bins)), bins_log)
axs[0].autoscale()

reg_widths = []
uhc_regs = []       # Ultra-high coverage regions
reg_start = -1
for idx, b in enumerate(bins):
    if b >= 3 * (bin_mean + bin_std):
        if reg_start == -1:
            reg_start = idx
    else:
        if reg_start > -1:
            reg_widths.append((idx - reg_start + 1) * nBinSize)
            uhc_regs.append({'start': reg_start * nBinSize, 'end': idx * nBinSize})
            axs[0].fill_between(range(reg_start, idx), bins_log[reg_start : idx], color='red')
            reg_start = -1
axs[0].axhline(math.log(bin_mean), color='green', linestyle='dotted')
axs[0].axhline(math.log(bin_mean + bin_std), color='orange', linestyle='dashed')

plt.sca(axs[0])
locs, labels = plt.yticks()
new_labels = []
for y in locs:
    new_labels.append(f'{int(math.exp(y)):,}')
locs = np.append(locs, math.log(bin_mean))
new_labels.append(f'Mean: {int(math.exp(locs[-1])):,}')
locs = np.append(locs, math.log(bin_mean + bin_std))
new_labels.append(f'Mean + 3SD: {int(math.exp(locs[-1])):,}')
plt.yticks(locs, new_labels, fontsize=16)

locs, labels = plt.xticks()
new_labels = []
for x in locs:
    new_labels.append(f'{x * nBinSize / 1_000_000:,}Mb')
new_labels[-1] = ''
plt.xticks(locs, new_labels, fontsize=16)
plt.xlim([0, len(bins)])
plt.ylim([0, max(bins_log)])

total_len = sum(reg_widths)
print(f'Total length of {len(reg_widths)} extremely highly covered regions: {total_len:,}bp ({total_len/len(dat) * 100:.2}%)')
print(f'Average size: {stats.mean(reg_widths):.2f}')

custom_lines = [Line2D([0], [0], color='green', linestyle='dotted', lw=4),
                Line2D([0], [0], color='orange', linestyle='dashed', lw=4),
                Line2D([0], [0], color='red', lw=4)]
axs[0].legend(custom_lines, ['Mean coverage', 'Mean + 3xSD', 'Coverage above mean+3SD'], fontsize=18)
plt.xlabel('Genomic position', fontsize=18)
plt.ylabel('Coverage', fontsize=18)
plt.title(f'Coverage of {strScaffoldID}', fontsize=20)

#plt.sca(axs[1])
#plt.xlabel('Genomic position', fontsize=18)
#plt.title(f'Coverage of {strScaffoldID}', fontsize=20)
#axs[1].set_frame_on(False)
#axs[1].tick_params(left=False)
#axs[1].tick_params(labelleft=False)

plt.savefig(strOutputPDF) 

repeats = []
nRead = 0
print(f'Loading the repeats from {strRepeatMasker}')
with(open(strRepeatMasker, 'r')) as hFile:
    for line in hFile.readlines():
        nRead += 1
        if line.startswith('#'):
            continue
        scafID, _, _, repeatStart, repeatEnd, rest = line.split('\t', 5)
        if scafID == strScaffoldID:
            repeats.append({'start': int(repeatStart), 'end': int(repeatEnd)})
print(len(repeats))
print(f'Read {nRead} lines')
tmp = []
for r in repeats:
    tmp.append([r['start'], r['end']])
repeats = pd.DataFrame(data=tmp, columns=['Start', 'End'])
print(repeats.shape)

# There are several possibilities
# A
# UHC region   |-------|       Repeat |------|           No overlap
#
# B
# UHC region   |---------------------|
# Repeat     |----| or  |---|  or |-----|   
#
# C
# UHC region |----| or  |---|  or |-----|
# Repeat        |---------------------|

nNonEmpty = 0
for r in uhc_regs:
    reg_len = r['end'] - r['start']
    tmp = repeats[repeats['Start'].between(r['start'], r['end'], inclusive=True) | 
                  repeats['End'].between(r['start'], r['end'], inclusive=True) |
                 (r['start'] >= repeats['Start']) & (r['end'] <= repeats['End'])]
    rep_len = 0
    nRep = 0
    if not tmp.empty:
        nRep = tmp.shape[0]
        for row in tmp.itertuples():
            rep_len += row.End - row.Start
    f = rep_len / reg_len * 100
    print(f"Region: {r['start']}-{r['end']} has {rep_len}bp ({f:.2f}%) in {nRep} repeats")

pickle.dump(uhc_regs, open(sys.argv[6], 'wb'))