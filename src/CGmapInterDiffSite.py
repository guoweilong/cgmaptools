#!/usr/bin/env python

"""
    cgmaptools - CGmapInterDiffSite.py

    Copyright (C) Weilong Guo
    Contact: Weilong Guo <guoweilong@126.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""

import sys
#import os
#import os.path
#import re

import gzip

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

try:
    from scipy.stats import chi2_contingency as chi2
    from scipy.stats import fisher_exact
except ImportError:
    sys.stderr.write("[Error] Cannot import \"scipy\" package. Have you installed it?")
    exit(-1)
#

# example 1
#chi2( [[19,24],[34,10]], correction=False )
#(9.9998158025027379, 0.0015655588405593997, 1, array([[ 26.1954023,  16.8045977],
#       [ 26.8045977,  17.1954023]]))
# example 2
#chi2( [[26,7],[36,2]], correction=True )
# example 3
#chi2( chi2([[180,14,120,65], [200,16,84,33], [380,30,204,98]]),correction=False )\

# By Fisher's Exact Test
# oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]]) >>> pvalue

def CGmapInterDiff (fn, MIN=0, MAX=100, method="chisq"):
    if fn is None:
        IN = sys.stdin
    else:
        try:
            if fn.endswith('.gz'):
                IN = gzip.open(fn, 'rt', encoding='UTF-8')
            else:
                IN = open(fn, 'r')
            #
        except IOError:
            print(f'\n[Error]\n\tFile cannot be open: {fn}')
            exit(-1)
        #
    #
    for line in IN:
        chr, nuc, pos, pattern, dinuc, methyl_1, NmC_1, NC_1, methyl_2, NmC_2, NC_2 = line.strip().split()
        NmC_1, NC_1, NmC_2, NC_2 = [int(NmC_1), int(NC_1), int(NmC_2), int(NC_2)]
        if (NC_1 < MIN or NC_2 < MIN) or (NC_1 > MAX or NC_2 > MAX) or ((NmC_1+NmC_2)==0):
            continue
        #
        if method == "chisq":
            pv = chi2([[NmC_1, NC_1], [NmC_2, NC_2]], correction=False)[1]
        elif method == "fisher":
            oddsratio, pv = fisher_exact([[NmC_1, NC_1-NmC_1], [NmC_2, NC_2-NmC_2]])
        else:
            pv = 1
        #
        print(f'{chr}\t{nuc}\t{pos}\t{pattern}\t{dinuc}\t{methyl_1}\t{methyl_2}\t{pv:.2e}')
    #
    if IN is not sys.stdin:
        IN.close()
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools dms [-i <CGmapInter>] [-m 5 -M 100] [-o output]\n" \
            "      (aka CGmapInterDiffSite)\n" \
            "Description: \n" \
            "  Get the differentially methylated sites for two samples.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2017-01-20\n" \
            "Input Format, same as the output of CGmapIntersect.py:\n" \
            "    Chr1  C  3541  CG  CG  0.8  4  5  0.4  4  10\n" \
            "Output Format:\n" \
            "   chr1	C	4654	CG	CG	0.92	1.00	8.40e-01\n" \
            "   chr1	C	4658	CHH	CC	0.50	0.00	3.68e-04\n" \
            "   chr1	G	8376	CG	CG	0.62	0.64	9.35e-01"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmapInter", help="File name for CGmapInter, STDIN if omitted", metavar="FILE")
    parser.add_option("-m", "--min", dest="min", help="min coverage [default : %default]", metavar="INT", default=0)
    parser.add_option("-M", "--max", dest="max", help="max coverage [default : %default]", metavar="INT", default=100)
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if omitted. Compressed output if end with .gz")
    parser.add_option("-t", "--test-method", dest="method", help="chisq, fisher [default : %default]", metavar="STRING",
                      default="chisq")
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz'):
            sys.stdout = gzip.open(options.outfile, 'wt', encoding='UTF-8')
        else:
            sys.stdout = open(options.outfile, 'w')
        #
    #
    if options.method not in ["chisq", "fisher"]:
        print(f'\n[Error]\n\tUnknown method : {options.method}')
        exit(-1)
    #
    CGmapInterDiff(options.CGmapInter, int(options.min), int(options.max), options.method)
#
# ===========================================
if __name__ == "__main__":
    main()
#
