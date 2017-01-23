#!/usr/bin/env python

"""
    cgmaptools - CGmapToRegion.py

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

import gzip


def Get_key (str) :
    # This function should be consistent with the one in Sort_chr_pos.py
    match = re.match(r"^chr(\d+)", str, re.I)
    if match :
        chr = int(match.group(1))
    else :
        chr = str
    return [chr, str]
    #re.match(r"^chr(\d+)", "chr019_random").group(1)

def CGmapToRegion (CGmap_fn, region_fn):
    try:
        if CGmap_fn :
            if CGmap_fn.endswith(".gz") :
                CGMAP = gzip.open(CGmap_fn, 'rb')
            else :
                CGMAP = open(CGmap_fn, 'r')
        else :
            CGMAP = sys.stdin
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", CGmap_fn
        exit(-1)
    #
    try:
        if region_fn.endswith(".gz") :
            REGION = gzip.open(region_fn, 'rb')
        else :
            REGION = open(region_fn, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", region_fn
        exit(-1)
    #
    line_c = CGMAP.readline()
    line_r = REGION.readline()
    #
    mC = 0
    count = 0
    NmC = 0
    NC = 0
    #
    try :
        chr_r, pos_r1, pos_r2 = line_r.strip().split()[0:3]
    except ValueError :
        print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % region_fn)
        exit(-1)
    #
    while line_c and line_r :
        #print "\t"+line_c, "\t"+line_r
        #print "========"
        try :
            chr_c, nuc_c, pos_c, pattern_c, dinuc_c, methyl_c, NmC_c, NC_c = line_c.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % CGmap_fn)
            exit(-1)
        try :
            chr_r, pos_r1, pos_r2 = line_r.strip().split()[0:3]
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % region_fn)
            exit(-1)
        #print chr_r, pos_r1, pos_r2
        key_c = Get_key(chr_c)
        key_r = Get_key(chr_r)
        if key_c < key_r :
            line_c = CGMAP.readline()
        elif key_c > key_r :
            if count > 0 :
                print "%s\t%s\t%s\t%.4f\t%d\t%.4f\t%d" % (chr_r, pos_r1, pos_r2, mC/count, count, float(NmC)/NC, NC)
            else :
                print "%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_r, pos_r1, pos_r2, count, NC)
            # init
            mC = 0
            count = 0
            NmC = 0
            NC = 0
            # do the next
            line_r = REGION.readline()
        else : # chr_c == chr_r
            if int(pos_c) < int(pos_r1) :
                line_c = CGMAP.readline()
            elif int(pos_c) > int(pos_r2) :
                if count > 0 :
                    print "%s\t%s\t%s\t%.4f\t%d\t%.4f\t%d" % (chr_r, pos_r1, pos_r2, mC/count, count, float(NmC)/NC, NC)
                else :
                    print "%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_r, pos_r1, pos_r2, count, NC)
                # init
                mC = 0
                count = 0
                NmC = 0
                NC = 0
                line_r = REGION.readline()
            else : # pos_r1 <= chr_c <= pos_r2
                NC_c = int(NC_c)
                if NC_c > 0 :
                    mC = mC + ( float(NmC_c) / NC_c )
                    count = count + 1
                    NmC = NmC + int(NmC_c)
                    NC = NC + int(NC_c)
                line_c = CGMAP.readline()
            #
        #
    # End for reading files
    #
    if count > 0 :
        print "%s\t%s\t%s\t%.2f\t%d\t%.2f\t%d" % (chr_r, pos_r1, pos_r2, mC/count, count, float(NmC)/NC, NC)
    else :
        print "%s\t%s\t%s\tNA\t%d\tNA\t%d" % (chr_r, pos_r1, pos_r2, count, NC)
    #
    REGION.close()
    if CGMAP is not sys.stdin:
        CGMAP.close()
    #
#
from optparse import OptionParser

# todo:
# * Provide the build-in sorting function

# ===========================================
def main():
    usage = "Usage: cgmaptools mtr [-i <CGmap>] -r <region> [-o <output>]\n" \
            "      (aka CGmapToRegion)\n" \
            "Description: Calculated the methylation levels in regions in two ways.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2017-01-20\n" \
            "Format of Region file:\n" \
            "  #chr    start_pos  end_pos\n" \
            "   chr1   8275       8429\n" \
            "Output file format:\n" \
            "   #chr  start_pos  end_pos  mean(mC)  #_C  #read(C)/#read(T+C)  #read(T+C)\n" \
            "   chr1   8275       8429     0.34     72         0.40             164\n" \
            "Note: The two input CGmap files should be sorted by Sort_chr_pos.py first.\n" \
            "      This script would not distinguish CG/CHG/CHH contexts."
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmapFile", help="File name end with .CGmap or .CGmap.gz. "
                                                   "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-r", dest="regionFile", help="Filename for region file, support *.gz", metavar="FILE")
    parser.add_option("-o", dest="outfile", default=None, help="To standard output if not specified.\n")
    #
    (options, args) = parser.parse_args()
    #
    if (options.regionFile is None) :
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
    CGmapToRegion(options.CGmapFile, options.regionFile)
#


# ===========================================
if __name__ == "__main__":
    main()
#
