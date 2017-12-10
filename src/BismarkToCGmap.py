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

# Guo, Weilong; guoweilong@126.com; 2017-08-18
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

def BismarkCGmap (Bismark_fn, CGmap_fn):
    # ================
    try:
        if Bismark_fn :
            if Bismark_fn.endswith(".gz") :
                BismarkF = gzip.open(Bismark_fn, "rb")
            else :
                BismarkF = open(Bismark_fn, 'r')
            #
        else :
            BismarkF = sys.stdin
        #
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", Bismark_fn
        exit(-1)
    #
    # ==================
    try:
        if CGmap_fn is not None:
            if CGmap_fn.endswith(".gz") :
                CGmapF = gzip.open(CGmap_fn, "wb")
            else :
                CGmapF = open(CGmap_fn, 'w')
            #
        else :
            CGmapF = sys.stdout
        #
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", CGmap_fn
        exit(-1)
    #
    cur_chr = ""
    #
    line = BismarkF.readline()
    #
    while line:
        try :
            chr, pos, strand, NmC, NnC, pattern, trinuc = line.strip().split()
        except ValueError :
            print( "\n[Error]:\n\t Your input file [ %s ] has wrong number of columns." % CGmap_fn)
            print( "\t You may check whether the input file is correct." )
            print( "\t The input shall be something like \"*.CpG_report.txt.gz\"." )
            exit(-1)
        #
        if strand == "+" :
            nuc = "C"
        else :
            nuc = "G"
        #
        dinuc = trinuc[0:2]
        AllC = int(NmC) + int(NnC)
        if AllC > 0 :
            ratio = float(NmC)/AllC
            CGmapF.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%d\n" % (chr, nuc, pos, pattern, dinuc, ratio, NmC, AllC))
        #
        line = BismarkF.readline()
        #
    #
    # End for reading files
    if CGmapF is not sys.stdout :
        CGmapF.close()
    #
    if BismarkF is not sys.stdin :
        BismarkF.close()
    #
#
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools convert cgmap2wig [-i <CGmap>] [-w <wig>] [-c <INT> -b <float>]\n" \
            "      (aka CGmapToWig)\n" \
            "Description: Generate WIG file from CGmap.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2017-08-18"
    #
    parser = OptionParser(usage)
    #
    parser.add_option("-i", dest="Bismark", default=None, help="The output file of Bismark"
                            "Input file name end with .gz is considered as compressed, "
                            "use STDIN when not specified", metavar="FILE")
    parser.add_option("-o", "--CGmap", dest="CGmap", default=None, help="Output in CGmap format, "
                            "use STDOUT if omitted (gzipped if end with \'.gz\')", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    BismarkCGmap(options.Bismark, options.CGmap)
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#
#
