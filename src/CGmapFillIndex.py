#!/usr/bin/env python

"""
    cgmaptools - CGmapFillIndex.py

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

#import math
#x=float('nan')
#math.isnan(x)

import gzip

def CGmapFillIndex(index, CGmap_list, tag_list, min_cov=1, max_cov=200) :
    try:
        if index :
            if index.endswith('.gz') :
                INDEX = gzip.open(index, 'rb')
            else :
                INDEX = open(index, 'r')
        else :
             INDEX = sys.stdin
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", index
        exit(-1)
    #
    CGmap_lst = CGmap_list.split(",")
    tag_lst = tag_list.split(",")
    if tag_lst != "" and len(CGmap_lst) != len(tag_lst) :
        print >> sys.stderr, "[Error]: lengths of two lists are not equal\n"
        exit(-1)
    #
    N = len(CGmap_lst)
    #
    # Initialization
    Site = dict()
    Site_lst = []
    for line in INDEX :
        # chr10   100005504
        pos = line.strip()
        Site[pos] = ["nan"]*N
        Site_lst.append(pos)
    #
    if index :
        INDEX.close()
    #
    # Read all the CGmap files
    for i in xrange(N) :
        fn = CGmap_lst[i]
        try :
            if fn.endswith('.gz') :
                INPUT = gzip.open(fn, 'rb')
            else :
                INPUT = open(fn, 'r')
        except IOError :
            print "\n[Error]:\n\t File cannot be open: ", CGmap_lst[i]
            exit(-1)
        #
        # chr1    C       4654    CG      CG      0.846153846154  11      13
        tokens = []
        for line in INPUT :
            tokens = line.strip().split()
            pos = "\t".join([tokens[0], tokens[2]])
            cov = int(tokens[7])
            if (pos in Site) and (min_cov <= cov <= max_cov) :
                Site[pos][i] = "%.2f" % float(tokens[5])
            #
        #
        INPUT.close()
    #
    # Write the output files
    if tag_lst != "" : # if not tag was specified, do not print the tag line
        print "\t".join( ["chr", "pos"] + tag_lst)
    N_site = len(Site_lst)
    for i in xrange(N_site) :
        pos = Site_lst[i]
        print "\t".join( [pos] + Site[pos] )
    #
#


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools mergelist tomatrix  [-i <index>] -f <IN1,IN2,..> -t <tag1,tag2,..> [-o output]\n" \
            "      (aka CGmapFillIndex)\n" \
            "Description: Fill methylation levels according to the Index file for CGmap files in list.\n" \
            "Contact: Guo, Weilong; guoweilong@126.com;\n" \
            "Last Updated: 2016-12-07\n" \
            "Index format Ex:\n" \
            "   chr10   100005504\n" \
            "Output format Ex:\n" \
            "   chr     pos     tag1    tag2    tag3\n" \
            "   Chr1    111403  0.30    nan     0.80\n" \
            "   Chr1    111406  0.66    0.40    0.60"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="index", help="TXT file, index file, use STDIN if omitted", metavar="FILE")
    parser.add_option("-f", dest="CGmap_list", help="List of (input) CGmap files (CGmap or CGmap.gz)", metavar="STRING")
    parser.add_option("-t", dest="tag_list",  metavar="STRING",
                      help="List of tags, same order with \'-f\' ", default="")
    parser.add_option("-c", dest="min_Cov", help="minimum coverage [default: %default]", metavar="INT", default=1)
    parser.add_option("-C", dest="max_Cov", help="maximum coverage [default: %default]", metavar="INT", default=200)
    parser.add_option("-o", dest="outfile", help="Output file name (gzipped if end with \'.gz\')", metavar="STRING")
    (options, args) = parser.parse_args()
    #
    if (options.index is None or options.CGmap_list is None or options.tag_list is None) :
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
    CGmapFillIndex(options.index, options.CGmap_list, options.tag_list, int(options.min_Cov), int(options.max_Cov))
#

# ===========================================
if __name__ == "__main__":
    main()

