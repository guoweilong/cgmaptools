#!/usr/bin/env python

"""
    cgmaptools - Sort_chr_pos.py

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
import re

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

import gzip
#import numpy as np
#import matplotlib.pyplot as plt

import datetime

def logm(message):
    log_message = "[%s] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print log_message
#

def SortMap (fn, chr_col=1, pos_col=2):
    chr_col_index = int(chr_col) - 1
    pos_col_index = int(pos_col) - 1
    try:
        if fn :
            if fn.endswith(".gz") :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
            #
        else :
            IN = sys.stdin
        #
    except IOError:
        sys.stderr.write("[Error] File cannot be open: %s\n" % fn)
        exit(-1)
    #
    lines = [ line.strip() for line in IN.readlines()]
    if fn :
        IN.close()
    #
    def Get_key (str) :
        # This function is linked to the one in CGmapToRegion.py
        tokens = str.split('\t')
        match = re.match(r"^chr(\d+)", tokens[chr_col_index], re.I)
        #chr = int(match.group(1)) if match else tokens[0]
        if match :
            chr = int(match.group(1))
        else :
            chr = tokens[chr_col_index]
        #
        try :
            pos = int(tokens[pos_col_index])
        except ValueError :
            sys.stderr.write("[Error] wrong content for position column.\n")
            exit(-1)
        #
        return [chr, tokens[chr_col_index], pos]
        #re.match(r"^chr(\d+)", "chr019_random").group(1)
    #
    lines.sort(key=lambda l: Get_key(l))
    #
    for i in lines :
        print("%s" % i)
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: %prog [-i <input>] [-c 1] [-p 3] [-o output]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2014-05-11\n" \
            "Last Update: 2016-12-07\n" \
            "Description: Sort the input files by chromosome and position.\n" \
            "     The order of chromosomes would be :\n" \
            "     \"chr1 chr2 ... chr11 chr11_random ... chr21 ... chrM chrX chrY\""
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="input",
                      help="File name end with .CGmap or .CGmap.gz. "
                           "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-c", "--chr", dest="chr_col",
                      help="The column of chromosome [default: %default]",
                      default=1, metavar="INT")
    parser.add_option("-p", "--pos", dest="pos_col",
                      help="The column of position [default: %default]",
                      default=2, metavar="INT")
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified")
    #
    (options, args) = parser.parse_args()
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    #
    #
    SortMap(options.input, options.chr_col, options.pos_col)
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
