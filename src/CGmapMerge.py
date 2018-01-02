#!/usr/bin/env python

"""
    cgmaptools - CGmapMerge.py

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

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

import gzip

def MergeCGmap (fn1, fn2):
    try:
        if fn1.endswith(".gz") :
            IN_1 = gzip.open(fn1, 'rb')
        else :
            IN_1 = open(fn1, 'r')
        #
    except IOError :
        print("\n[Error]:\n\t File cannot be open: %s" % fn1)
        exit(-1)
    #
    try:
        if fn2 :
            if fn2.endswith(".gz") :
                IN_2 = gzip.open(fn2, 'rb')
            else :
                IN_2 = open(fn2, 'r')
            #
        else :
            IN_2 = sys.stdin
        #
    except IOError:
        print("\n[Error]:\n\t File cannot be open: %s" % fn2 )
        exit(-1)
    #
    line_1 = IN_1.readline()
    line_2 = IN_2.readline()
    chr_pre = ""
    #
    while line_1 and line_2 :
        chr_1, nuc_1, pos_1, pattern_1, dinuc_1, methyl_1, NmC_1, NC_1 = line_1.strip().split()
        chr_2, nuc_2, pos_2, pattern_2, dinuc_2, methyl_2, NmC_2, NC_2 = line_2.strip().split()
        if chr_1 != chr_2 :
            if chr_2 != chr_pre:
                print line_1.strip()
                line_1 = IN_1.readline()
            elif chr_1 != chr_pre :
                print line_2.strip()
                line_2 = IN_2.readline()
            #
        else : # chr_1 == chr_2
            chr_pre = chr_1
            if int(pos_1) < int(pos_2) :
                print line_1.strip()
                line_1 = IN_1.readline()
            elif int(pos_1) > int(pos_2) :
                print line_2.strip()
                line_2 = IN_2.readline()
            else :
                if nuc_1 != nuc_2 or pattern_1 != pattern_2 or dinuc_1 != dinuc_2 :
                    sys.stderr.write("Warning: Inconsistent information:")
                    sys.stderr.write("%s | %s" % (line_1, line_2) )
                #
                NmC = int(NmC_1) + int(NmC_2)
                NC = int(NC_1) + int(NC_2)
                methyl = float(NmC)/NC
                print("\t".join([chr_1, nuc_1, pos_1, pattern_1, dinuc_1, "%.2f" % methyl, "%d" % NmC, "%d" % NC]) )
                line_1 = IN_1.readline()
                line_2 = IN_2.readline()
            #
        #
    #
    # End for reading files
    while line_1 :
        print line_1.strip()
        line_1 = IN_1.readline()
    #
    while line_2 :
        print line_2.strip()
        line_2 = IN_2.readline()
    #
    IN_1.close()
    if IN_2 is not sys.stdin:
        IN_2.close()
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools merge2 cgmap -1 <CGmap_1> -2 <CGmap_2> [-o <output>]\n" \
            "      (aka CGmapMerge)\n" \
            "Description: Merge two CGmap files together.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-01-02\n" \
            "Note: The two input CGmap files should be sorted in the same order first.\n"
    parser = OptionParser(usage)
    parser.add_option("-1", dest="CGmap_1",
                      help="File name end with .CGmap or .CGmap.gz", metavar="FILE")
    parser.add_option("-2", dest="CGmap_2",
                      help="If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-o", dest="outfile", default=None,
                      help="CGmap, output file. Use STDOUT if omitted "
                           "(gzipped if end with \'.gz\').")
    #
    (options, args) = parser.parse_args()
    #
    if (options.CGmap_1 is None) :
        parser.print_help()
        exit(-1)
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    MergeCGmap(options.CGmap_1, options.CGmap_2)
#


# ===========================================
if __name__ == "__main__":
    main()
#
