#!/usr/bin/env python

"""
    cgmaptools - CGmapToWig.py

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

# Guo, Weilong; guoweilong@126.com; 2013-12-11
import sys
#import os
#import os.path
#import re

""" CGmap file
chr1    G   3000851 CHH CC  0.1 1   10
chr1    C   3001624 CHG CA  0.0 0   9
chr1    C   3001631 CG  CG  1.0 5   5
"""

""" WIG file
variableStep chrom=chr1
3000419 0.000000
3000423 -0.2
3000440 0.000000
3000588 0.5
3000593 -0.000000
"""


import gzip

def CGmapToWig (CGmap_fn, WIG_fn, coverage=1, base=0):
    base = float(base)
    coverage = int(coverage)
    try:
        if CGmap_fn :
            if CGmap_fn.endswith(".gz") :
                CGmap = gzip.open(CGmap_fn, "rb")
            else :
                CGmap = open(CGmap_fn, 'r')
            #
        else :
            CGmap = sys.stdin
        #
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", CGmap_fn
        exit(-1)
    #
    try:
        if WIG_fn is not None:
            if WIG_fn.endswith(".gz") :
                WIG = gzip.open(WIG_fn, "wb")
            else :
                WIG = open(WIG_fn, 'w')
            #
        else :
            WIG = sys.stdout
        #
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", WIG_fn
        exit(-1)
    #
    cur_chr = ""
    line = CGmap.readline()
    while line:
        try :
            chr, nuc, pos, pattern, dinuc, methyl, un_C, all_C = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % CGmap_fn)
            exit(-1)
        #
        if cur_chr != chr :
            cur_chr = chr
            WIG.write("variableStep chrom=" + chr + "\n")
        if methyl != "na" and int(all_C) >= coverage :
            methyl=float(methyl)+base
            if nuc == "C" :
                WIG.write("%s\t%.2f\n" % (pos, methyl) )
            elif nuc == "G":
                WIG.write("%s\t-%.2f\n" % (pos, methyl) ) # Negative value for Minus strand
            #
        #
        line = CGmap.readline()
    #
    # End for reading files
    if WIG is not sys.stdout  :
        WIG.close()
    #
    if CGmap is not sys.stdin:
        CGmap.close()
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools convert cgmap2wig [-i <CGmap>] [-w <wig>] [-c <INT> -b <float>]\n" \
            "      (aka CGmapToWig)\n" \
            "Description: Generate WIG file from CGmap.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2016-12-07"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap", default=None,
                      help="Input file name end with .CGmap or .CGmap.gz, "
                           "use STDIN when not specified", metavar="FILE")
    parser.add_option("-w", dest="WIG", default=None,
                      help="Output, use STDOUT if omitted (gzipped if end with \'.gz\')",
                      metavar="FILE")
    parser.add_option("-c", dest="coverage", default=1,
                      help="The minimum coverage of the sites [default: %default]",
                      metavar="INT")
    parser.add_option("-b", dest="base", default=0,
                      help="The base for adding to distinguish '0' and 'Nan' [default: %default]",
                      metavar="FLOAT")
    #
    (options, args) = parser.parse_args()
    #
    CGmapToWig(options.CGmap, options.WIG, options.coverage, options.base )
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#
