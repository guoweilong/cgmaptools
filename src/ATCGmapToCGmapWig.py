#!/usr/bin/env python

"""
    cgmaptools - ATCGmapToCGmapWig.py

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

""" ATCGmap File
chr1    T   3009410 --  --  0   10  0   0   0   0   0   0   0   0   na
chr1    C   3009411 CHH CC  0   10  0   0   0   0   0   0   0   0   0.0
chr1    C   3009412 CHG CC  0   10  0   0   0   0   0   0   0   0   0.0
chr1    C   3009413 CG  CG  0   10  50  0   0   0   0   0   0   0   0.83
"""

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

def ATCGmapToCGmapWig (ATCGmap_fn, CGmap_fn, WIG_fn):
    try:
        if ATCGmap_fn :
            if ATCGmap_fn.endswith(".gz") :
                ATCGmap = gzip.open(ATCGmap_fn, "rb")
            else :
                ATCGmap = open(ATCGmap_fn, 'r')
            #
        else :
            ATCGmap = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % ATCGmap_fn)
        exit(-1)
    #
    if (CGmap_fn is not None) :
        if CGmap_fn.endswith('.gz') :
            sys.stdout = gzip.open(CGmap_fn, 'wb')
        else :
            sys.stdout = open(CGmap_fn, 'w')
        #
    #
    if (WIG_fn is not None):
        try:
            if WIG_fn.endswith(".gz") :
                WIG = gzip.open(WIG_fn, "wb")
            else :
                WIG = open(WIG_fn, 'w')
            #
        except IOError:
            print ("\n[Error]:\n\t File cannot be open: %s" % WIG_fn)
            exit(-1)
        #
    #
    cur_chr = ""
    #
    line = ATCGmap.readline()
    #
    while line:
        try :
            chr, nuc, pos, pattern, dinuc, WA, WT, WC, WG, WN, CA, CT, CC, CG, CN, methyl = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % ATCGmap_fn)
            exit(-1)
        #
        if cur_chr != chr :
            cur_chr = chr
            if WIG_fn is not None :
                WIG.write("variableStep chrom=" + chr + "\n")
            #
        #
        if methyl != "na" :
            if nuc == "C" :
                print("\t".join([chr, nuc, pos, pattern, dinuc, methyl, WC, "%d" % (int(WT)+int(WC)) ]))
                if WIG_fn :
                    WIG.write(pos + " " + methyl + "\n")
                #
            elif nuc == "G":
                print("\t".join([chr, nuc, pos, pattern, dinuc, methyl, CG, "%d" % (int(CG)+int(CA)) ]))
                if WIG_fn :
                    WIG.write(pos + " -" + methyl + "\n") # Negative value for Minus strand
                #
            #
        #
        line = ATCGmap.readline()
        #
    #
    # End for reading files
    #
    if WIG_fn is not None :
        WIG.close()
    #
    if ATCGmap is not sys.stdin:
        ATCGmap.close()
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools convert atcgmap2cgmap [-i <ATCGmap>] [-c <CGmap>] [-w <wig>]\n" \
            "       (aka ATCGmapToCGmapWig)\n" \
            "Description: Generate CGmap or WIG from ATCGmap.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-01-02"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="ATCGmap", default=None,
                      help= "Input file name end with .ATCGmap or .ATCGmap.gz, "
                            "use STDIN when not specified", metavar="FILE")
    parser.add_option("-c", dest="CGmap", default=None,
                      help= "Output, gzip if end with \".gz\". "
                            "use STOUT if not specified",  metavar="FILE")
    parser.add_option("-w", dest="WIG", default=None, help="Output, gzip if end with \".gz\"",
                      metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    ATCGmapToCGmapWig(options.ATCGmap, options.CGmap, options.WIG)
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#
