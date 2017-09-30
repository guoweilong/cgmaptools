#!/usr/bin/env python

"""
    cgmaptools - CGmapSlitByChr.py

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

# Start from 2015-05-22

import sys
#import os
#import os.path
#import re

""" CGmap file
chr1    G   3000851 CHH CC  0.1 1   10
chr1    C   3001624 CHG CA  0.0 0   9
chr1    C   3001631 CG  CG  1.0 5   5
"""

#
import gzip
#
def CGmapSplitByChr (CGmap_fn, output_prefix, output_suffix):
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
        sys.stderr.write( "[Error] file %s cannot be open." % CGmap_fn)
        exit(-1)
    #
    line = CGmap.readline()
    tokens = line.strip().split()
    chr = tokens[0]
    if output_suffix.endswith(".gz") :
        OUT = gzip.open (output_prefix + "." + chr + "." + output_suffix, "wb")
    else :
        OUT = open (output_prefix + "." + chr + "." + output_suffix, "w")
    #
    cur_chr = chr
    while line:
        tokens = line.strip().split()
        chr = tokens[0]
        if chr != cur_chr :
            OUT.close()
            cur_chr = chr
            if output_suffix.endswith(".gz") :
                OUT = gzip.open (output_prefix + "." + chr + "." + output_suffix, "wb")
            else :
                OUT = open (output_prefix + "." + chr + "." + output_suffix, "w")
            #
            OUT.write(line)
        else :
            OUT.write(line)
        #
        line = CGmap.readline()
        #
    #
    if CGmap is not sys.stdin:
        CGmap.close()
    #
    OUT.close()
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools split -i <input> -p <prefix[.chr.]> -s <[.chr.]suffix>\n" \
            "      (aka CGmapSplitByChr)\n" \
            "Description: Split the files by each chromosomes. \n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2016-12-07"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap", default=None,
                      help="Input file, CGmap or ATCGmap foramt, "
                           "use STDIN when not specified."
                           "(gzipped if end with \'gz\').", metavar="FILE")
    parser.add_option("-p", dest="output_prefix", default=0,
                      help="The prefix for output file", metavar="STRING")
    parser.add_option("-s", dest="output_suffix", default=0,
                      help="The suffix for output file "
                           "(gzipped if end with \'gz\').", metavar="STRING")
    #
    (options, args) = parser.parse_args()
    #
    CGmapSplitByChr(options.CGmap, options.output_prefix, options.output_suffix )
#
# ===========================================
if __name__ == "__main__":
    main()
#
