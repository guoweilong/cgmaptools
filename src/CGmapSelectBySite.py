#!/usr/bin/env python

"""
    cgmaptools - CGmapSelectBySite.py

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

def CGmapInIndex(index, CGmap, reverse) :
    try:
        if index.endswith(".gz") :
            INDEX = gzip.open(index, 'rb')
        else :
            INDEX = open(index, 'r')
        #
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", index
        exit(-1)
    #
    # Initialization
    Site_lst = []
    for line in INDEX :
        # chr10   100005504
        pos = line.strip()
        Site_lst.append(pos)
    #
    INDEX.close()
    #
    # chr1    C       4654    CG      CG      0.846153846154  11      13
    try :
        with (gzip.open(CGmap, 'rb') if CGmap.endswith(".gz") else open(CGmap, 'r')  )  if CGmap else sys.stdin as IN:
            for line in IN :
                line = line.strip()
                tokens = line.split()
                pos = "\t".join([tokens[0], tokens[2]])
                if not reverse :
                    if pos in Site_lst :
                        print line
                    #
                else:
                    if pos not in Site_lst :
                        print line
                    #
                #
            #
        #
        if CGmap :
            IN.close()
        #
    except IOError :
        print "\n[Error]:\n\t File cannot be open: ", CGmap
        exit(-1)
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools select site -i <index> [-f <CGmap/ATCGmap>] [-r] [-o output]\n" \
            "      (aka CGmapSelectBySite)\n" \
            "Description: Select lines from input CGmap/ATCGmap in index or reverse.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2016-12-07\n" \
            "Index format example:\n" \
            "   chr10   100504\n" \
            "   chr10   103664"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="index",
                      help= "Name of Index file required "
                            "(gzipped if end with \'.gz\').", metavar="FILE")
    parser.add_option("-r", action="store_true", dest="reverse",
                      help= "reverse selected, remove site in index if specified", default = False)
    parser.add_option("-f", dest="infile",
                      help= "Input CGmap/ATCGmap files. Use STDIN if not specified", metavar="STRING")
    parser.add_option("-o", dest="outfile",
                      help= "CGmap, Output file name (gzipped if end with \'.gz\').", metavar="STRING")
    (options, args) = parser.parse_args()
    #
    if (options.index is None) :
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
    CGmapInIndex(options.index, options.infile, options.reverse)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
