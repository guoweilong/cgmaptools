#!/usr/bin/env python

"""
    cgmaptools - FragRegFromBED.py

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

#import os
import gzip
import sys
#import re

def error(msg):
    sys.stderr.write(f'ERROR: {msg}')
    exit(1)
#
from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools bed2fragreg [-i <BED>] [-n <N>] [-F <50,50,..> -T <50,..>] [-o output]\n" \
            "      (aka FragRegFromBED)\n" \
            "Description: Generate fragmented regions from BED file.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-05-02\n" \
            "   Split input region into N bins, get fragments from 5' end and 3' end.\n" \
            "Input Ex:\n" \
            "   chr1   1000    2000   +\n"\
            "   chr2   9000    8000   -\n"\
            "Output Ex:\n" \
            "   chr1   +   940  950  1000 1200 1400 1600 1800 1850\n" \
            "   chr2   -   9060 9050 9000 8800 8600 8400 8200 8150\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile", help="BED format, STDIN if omitted", metavar="FILE")
    parser.add_option("-F", dest="FiveMerEnd",
                      help="List of region lengths in upstream of 5\' end, Ex: 10,50. "
                           "List is from 5\'end->3\'end",
                      metavar="INT_list", default="")
    parser.add_option("-T", dest="ThreeMerEnd",
                      help="List of region lengths in downstream of 3\' end, Ex: 40,20. "
                           "List is from 5\'end->3\'end",
                      metavar="INT_list", default="")
    parser.add_option("-n", dest="nbins", help="Number of bins to be equally split [Default:%default]\n",
                      metavar="INT", default=1)
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if omitted. Compressed output if end with .gz")
    (options, args) = parser.parse_args()
    #
    if (options.infile is None):
        INPUT = sys.stdin
    else:
        try:
            if options.infile.endswith('.gz'):
                INPUT =  gzip.open(options.infile, 'rt', encoding='UTF-8')
            else:
                INPUT =  open(options.infile, 'r')
            #
        except IOError:
            print(f'[Error] Cannot find input file : {options.infile} !')
            exit(-1)
        #
    #
    if (options.outfile is not None):
        if options.outfile.endswith('.gz'):
            sys.stdout = gzip.open(options.outfile, 'wt', encoding='UTF-8')
        else:
            sys.stdout = open(options.outfile, 'w')
        #
    #
    #
    # works for python2 and python3
    try:
        xrange
    except NameError:
        xrange = range
    #
    # FiveMerEND
    if options.FiveMerEnd == "":
        FiveMerLst=[]
    else:
        FiveMerLst=[int(i) for i in options.FiveMerEnd.split(",")]
    #
    FiveMerPos = []
    pos_tmp = 0
    for i in FiveMerLst[::-1]:
        pos_tmp += i
        FiveMerPos = [pos_tmp] + FiveMerPos
    #
    # ThreeMerEND
    if options.ThreeMerEnd == "":
        ThreeMerLst=[]
    else:
        ThreeMerLst=[int(i) for i in options.ThreeMerEnd.split(",")]
    #
    ThreeMerPos = []
    pos_tmp = 0
    for i in ThreeMerLst:
        pos_tmp += i
        ThreeMerPos = ThreeMerPos + [pos_tmp]
    #
    nbins = int(options.nbins)
    # print the header line
    print("\t".join(["#chr", "chr"] + ["U%d" %(i+1) for i in range(len(FiveMerLst))]
           + ["R%d" % (i + 1) for i in range(nbins)]
           + ["D%d" % (i + 1) for i in range(len(ThreeMerLst))] + ['End']))
    # Example for BED line
    # chr6    -       33661860        33669093
    for line in INPUT:
        #print line.strip().split()
        [chr, Left, Right, strand]=line.strip().split()
        Left = int(Left)
        Right = int(Right)
        if strand == '+':
            body_regions = [(Left+(Right-Left)*i/nbins) for i in xrange(nbins+1)]
            LeftEnd_regions = [(Left-i) for i in FiveMerPos]
            RightEnd_regions = [(Right+i) for i in ThreeMerPos]
            print("\t".join([chr, strand]) + "\t" + "\t".join(str(int((i+abs(i))/2)) for i in (LeftEnd_regions + body_regions + RightEnd_regions)))
            #
        else:
            body_regions = [(Left+(Right-Left)*i/nbins) for i in xrange(nbins+1)]
            LeftEnd_regions = [(Left-i) for i in ThreeMerPos]
            RightEnd_regions = [(Right+i) for i in FiveMerPos]
            print("\t".join([chr, strand]) + "\t" + "\t".join(str(int((i+abs(i))/2)) for i in (RightEnd_regions + body_regions[::-1] + LeftEnd_regions)))
            #
        #
    #
    if options.infile is not None:
        INPUT.close()
    #
#
# ===========================================
if __name__ == "__main__":
    main()
#
