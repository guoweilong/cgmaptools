#!/usr/bin/env python

"""
    cgmaptools - MergeListOfCGmap.py

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
chr1    T   3009410 --  --  0   10  0   0   0   0   0   0   0   0   na
chr1    C   3009411 CHH CC  0   10  0   0   0   0   0   0   0   0   0.0
chr1    C   3009412 CHG CC  0   10  0   0   0   0   0   0   0   0   0.0
chr1    C   3009413 CG  CG  0   10  50  0   0   0   0   0   0   0   0.83
"""



import gzip

def MergeListOfATCGmap (fn_lst):
    ATCGmap_lst = fn_lst.split(",")
    #
    Methylome = {}
    for fn in ATCGmap_lst :
        try :
            if fn.endswith('.gz') :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
            #
        except IOError :
            print("\n[Error]:\n\t File cannot be open: %s" % fn )
            exit(-1)
        #
        for line in IN :
            try :
                chr, nuc, pos, pattern, dinuc, WA, WT, WC, WG, WN, CA, CT, CC, CG, CN, methyl = line.strip().split()
            except ValueError :
                print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
                exit(-1)
            #
            index = chr + "\t" + pos
            [WA, WT, WC, WG, WN, CA, CT, CC, CG, CN] = [int(WA), int(WT), int(WC), int(WG), int(WN), int(CA), int(CT), int(CC), int(CG), int(CN)]
            if index not in Methylome :
                Methylome[index] = [ "\t".join([nuc, pattern, dinuc]), [WA, WT, WC, WG, WN, CA, CT, CC, CG, CN] ]
            else :
                Methylome[index][1] = [i+j for i,j in zip(Methylome[index][1], [int(WA), int(WT), int(WC), int(WG), int(WN), int(CA), int(CT), int(CC), int(CG), int(CN)]) ]
            #
        #
        IN.close()
    #
    def Get_key (str) :
        # This function is linked to the one in CGmapToRegion.py
        tokens = str.split()
        match = re.match(r"^chr(\d+)", tokens[0], re.I)
        #chr = int(match.group(1)) if match else tokens[0]
        if match :
            chr = int(match.group(1))
        else :
            chr = tokens[0]
        #
        pos = int(tokens[1])
        return [chr, tokens[0], pos]
        #re.match(r"^chr(\d+)", "chr019_random").group(1)
    #
    for index in sorted(Methylome.keys(), key=lambda l: Get_key(l)) :
        [chr, pos] = index.split("\t")
        [nuc, pattern, dinuc] = Methylome[index][0].split()
        [WA, WT, WC, WG, WN, CA, CT, CC, CG, CN] = Methylome[index][1]
        methyl = 'na'
        if nuc == "C" :
            if WT+WC>0 :
                methyl = "%.2f" % (float(WC)/(WT+WC))
            #
        elif nuc == "G" :
            if CA+CG>0 :
                methyl = "%.2f" % (float(CG)/(CA+CG))
            #
        #
        print("%s\t%c\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s" %
              (chr, nuc, pos, pattern, dinuc, WA, WT, WC, WG, WN, CA, CT, CC, CG, CN, methyl) )
        #
    #
#

def MergeListOfCGmap (fn_lst):
    CGmap_lst = fn_lst.split(",")
    #
    Methylome = {}
    for fn in CGmap_lst :
        try :
            if fn.endswith('.gz') :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
            #
        except IOError :
            print("\n[Error]:\n\t File cannot be open: %s" % fn )
            exit(-1)
        #
        for line in IN :
            try :
                chr, nuc, pos, pattern, dinuc, methyl, MC, NC = line.strip().split()
            except ValueError :
                print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
                exit(-1)
            #
            index = chr + "\t" + pos
            [MC, NC] = [int(MC), int(NC)]
            if index not in Methylome :
                Methylome[index] = [ "\t".join([nuc, pattern, dinuc]), [MC, NC] ]
            else :
                Methylome[index][1] = [i+j for i,j in zip(Methylome[index][1], [int(MC), int(NC)]) ]
            #
        #
        IN.close()
    #
    def Get_key (str) :
        # This function is linked to the one in CGmapToRegion.py
        tokens = str.split()
        match = re.match(r"^chr(\d+)", tokens[0], re.I)
        if match :
            chr = int(match.group(1))
        else :
            chr = tokens[0]
        #
        pos = int(tokens[1])
        return [chr, tokens[0], pos]
    #
    for index in sorted(Methylome.keys(), key=lambda l: Get_key(l)) :
        [chr, pos] = index.split("\t")
        [nuc, pattern, dinuc] = Methylome[index][0].split()
        [MC, NC] = Methylome[index][1]
        methyl = 'na'
        if NC>0 :
            methyl = "%.2f" % (float(MC)/NC)
        #
        print("%s\t%c\t%s\t%s\t%s\t%s\t%d\t%d" %  (chr, nuc, pos, pattern, dinuc, methyl, MC, NC) )
        #
    #
#


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools mergelist tosingle -i f1,f2,..,fn [-o <output>]\n" \
            "      (aka MergeListOfCGmap)\n" \
            "Description: Merge multiple CGmap/ATCGmap files into one.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2016-12-07\n" \
            "Note: Large memory is needed.\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile_lst",
                      help="List of input files; gzipped file ends with '.gz'", metavar="FILE")
    parser.add_option("-f", dest="format",
                      help="cgmap or atcgmap [Default: %default]", default = "cgmap", metavar="FILE")
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified; gzipped file if end with '.gz'")
    #
    (options, args) = parser.parse_args()
    #
    if (options.infile_lst is None) :
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
    if options.format == "cgmap" :
        MergeListOfCGmap(options.infile_lst)
    else :
        MergeListOfATCGmap(options.infile_lst)
    #
#
# todo : automatic detect CGmap format or ATCGmap format


# ===========================================
if __name__ == "__main__":
    main()
#
