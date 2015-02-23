#!/usr/bin/python

import sys, os
import string
import random

'''
Take a bed file of boundaries (in this case, compartments)
and generates a series of bins around them, for use with
bigWigAverageOverBed
'''

try:
    fname = sys.argv[1]
    dist  = int(sys.argv[2])
except:
    print "\nCall:\n\tpython binAroundBed.py <input bed> <+/-dist> \n" +\
        "\te.g. python binAroundBed.py h1_tbounds.bed 500000\n"
    sys.exit()
# i.e. h1_superBounds.bed | h1_tbounds.bed | h1_compart_bounds.bed

if dist == 1500000:
    # compartment bounds
    res   = 100000
    steps = 30
elif dist == 500000:
    # TAD bounds
    res   = 40000
    steps = 25
else:
    sys.exit("Distance (<1>) should be 1.5 Mb or 100 kb")
    

csize = {}
with open(os.path.expanduser("../data/text/hg19.chrom.sizes.txt"), "r") as cs:
    for l in cs:
        c = l.split()
        csize[c[0]] = int(c[1])

with open(fname, "r") as f:
    for line in f:
        cols = line.split()
        csome = cols[0]
        s = int(cols[1])
        e = int(cols[2])
        mid = (e + s) / 2
        name = cols[3]
        # +/- 500kb, 50 * 20kb bins OR 100 * 10kb bins OR 25 * 40?
        start = mid - dist
        end = mid + dist
        if end > (csize[csome] - res):
            pass
        elif not start < 0:
            for b in range(1,(steps + 1)):
                snum = start + res * b
                print "\t".join([csome, str(snum), str(snum+res),
                            csome + "_" + str(snum) + "_" \
                            + "".join(random.sample(string.ascii_letters, 7))])

# then run something like:
# > bedtools intersect -a h1_sbounds_bins.bed
#     -b ctcf/wgEncodeOpenChromChipH1hescCtcfPk.narrowPeak -c > s1bounds_isect.out
#
