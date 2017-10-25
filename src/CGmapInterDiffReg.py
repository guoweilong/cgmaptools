#!/usr/bin/env python

"""
    cgmaptools - CGmapInterDiffReg.py

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

try :
    #from scipy.stats import ttest_1samp
    from scipy.stats import ttest_ind
except ImportError :
    sys.stderr.write("[Error] Cannot import \"scipy\" package. Have you installed it?")
    exit(-1)
#

try :
    #from numpy import mean
    from numpy import seterr
except ImportError :
    sys.stderr.write("[Error] Cannot import \"numpy\" package. Have you installed it?")
    exit(-1)
#

seterr(divide='ignore', invalid='ignore')
#from scipy.stats import chi2_contingency as chi2
#from scipy.stats import fisher_exact
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


from math import fsum
def average(x):
    return fsum(x)/float(len(x)) if x else 0
#

def CGmapInterDiffRegion (fn, minCov=0, maxCov=100, minStep=100, maxStep=500, minNSite = 2):
    if fn is None :
        IN = sys.stdin
    else :
        try:
            if fn.endswith('.gz') :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
        except IOError:
            print "\n[Error]\n\tFile cannot be open: ", fn
            exit(-1)
    #
    pre_pos = 0
    start_pos = 0
    lst_1 = []
    lst_2 = []
    [Tstat, PV] = [0, 1]
    pre_chr = ""
    #
    for line in IN:
        if len(line) == 0 :
            continue
        #
        try:
            chr, nuc, pos, pattern, dinuc, methyl_1, NmC_1, NC_1, methyl_2, NmC_2, NC_2 = line.strip().split()
        except ValueError:
            sys.stderr.write("An invalid line is skipped: %s" % line)
            continue
        #
        NmC_1, NC_1, NmC_2, NC_2 = [int(NmC_1), int(NC_1), int(NmC_2), int(NC_2)]
        # Coverage not valid
        if (NC_1 < minCov or NC_2 < minCov) or (NC_1 > maxCov or NC_2 > maxCov) or ((NmC_1+NmC_2)==0) :
            continue
        # Check position
        pos = int(pos)
        if abs(pos-pre_pos) > minStep or abs(pos - start_pos) > maxStep:
            # check the old fragment
            N_site = len(lst_1)
            if N_site >= minNSite and (start_pos<pre_pos):
                #[Tstat, PV] = ttest_1samp( [i-j for i, j in zip(lst_1, lst_2) ], 0 )
                mean_1 = average(lst_1)
                mean_2 = average(lst_2)
                [Tstat, PV] = ttest_ind( lst_1, lst_2 )
                print( "%s\t%d\t%d\t%.4f\t%.2e\t%.4f\t%.4f\t%d"
                       % (pre_chr, start_pos, pre_pos, Tstat, PV, mean_1, mean_2, N_site) )
                #
            # start a new fragment
            lst_1 = [ float(methyl_1) ]
            lst_2 = [ float(methyl_2) ]
            start_pos = pos
        else :
            # for the new fragment
            lst_1.append( float(methyl_1) )
            lst_2.append( float(methyl_2) )
        #
        pre_chr = chr
        pre_pos = pos
    # End-of-for
    if len(lst_1) >=2 :
        #[Tstat, PV] = ttest_1samp( [i-j for i, j in zip(lst_1, lst_2) ], [0])
        mean_1 = average(lst_1)
        mean_2 = average(lst_2)
        [Tstat, PV] = ttest_ind( lst_1, lst_2 )
        N_site = len(lst_1)
        if N_site >= minNSite and (start_pos < pre_pos) :
            print( "%s\t%d\t%d\t%.4f\t%.2e\t%.4f\t%.4f\t%d"
                % (pre_chr, start_pos, pre_pos, Tstat, PV, mean_1, mean_2, N_site) )
            #
        #
    #
    if IN is not sys.stdin:
        IN.close()
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools dmr [-i <CGmapInter>] [-m 5 -M 100] [-o output]\n" \
            "      (aka CGmapInterDiffReg)\n" \
            "Description: \n" \
            "  Get the differentially methylated regions using dynamic fragment strategy.\n" \
            "Author:  Guo, Weilong; guoweilong@126.com; \n" \
            "Last Updated: 2017-09-28\n" \
            "Input Format, same as the output of CGmapIntersect.py:\n" \
            "   chr1  C  3541  CG  CG  0.8  4  5  0.4  4  10\n" \
            "Output Format, Ex:\n" \
            "  #chr	    start   end	    t       pv          mC_A    mC_B    N_site\n" \
            "   chr1    1004572 1004574 inf     0.00e+00    0.1100	0.0000  20\n" \
            "   chr1    1009552 1009566 -0.2774 8.08e-01    0.0200	0.0300  15\n" \
            "   chr1    1063405 1063498 0.1435  8.93e-01    0.6333	0.5733  5\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmapInter", help="File name for CGmapInter, STDIN if omitted", metavar="FILE")
    parser.add_option("-c", "--minCov", dest="minCov", help="min coverage [default : %default]", metavar="INT", default=4)
    parser.add_option("-C", "--maxCov", dest="maxCov", help="max coverage [default : %default]", metavar="INT", default=500)
    parser.add_option("-s", "--minStep", dest="minStep", help="min step in bp [default : %default]", metavar="INT", default=100)
    parser.add_option("-S", "--maxStep", dest="maxStep", help="max step in bp [default : %default]", metavar="INT", default=1000)
    parser.add_option("-n", "--minNSite", dest="minNSite", help="min N sites [default : %default]", metavar="INT", default=5)
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if omitted. Compressed output if end with .gz")
    #parser.add_option("-t", "--test-method", dest="method", help="chisq, fisher [default : %default]", metavar="STRING",
    #                  default="chisq")
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    minCov = int(options.minCov)
    maxCov = int(options.maxCov)
    minStep = int(options.minStep)
    maxStep = int(options.maxStep)
    minNSite = int(options.minNSite)
    #
    CGmapInterDiffRegion (options.CGmapInter, minCov, maxCov, minStep, maxStep, minNSite)
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#
