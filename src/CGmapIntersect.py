#!/usr/bin/env python

"""
    cgmaptools - CGmapIntersect.py

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

import gzip

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

def CheckCtx(context, dinuc, option_ctx) :
    if option_ctx == '' :
        return True
    if option_ctx in ['CA', 'CT', 'CG', 'CC'] :
        if dinuc == option_ctx :
            return True
        else :
            return False
        #
    if option_ctx in ['CHG', 'CHH'] :
        if context == option_ctx :
            return True
        else :
            return False
        #
    if option_ctx == 'CH' :
        if context in ['CHG', 'CHH'] :
            return True
        else :
            return False
        #
    if option_ctx == 'CW' :
        if context in ['CA', 'CT'] :
            return True
        else :
            return False
        #
    return False
#

def CGmapIntersect (fn1, fn2, ctx = ""):
    try:
        if fn1.endswith('.gz') :
            IN_1 = gzip.open(fn1, 'rb')
        else :
            IN_1 = open(fn1, 'r')
    except IOError:
        sys.stderr.write("[Error] File %s cannot be open." % fn1)
        exit(-1)
    #
    try:
        if fn2 :
            if fn2.endswith('.gz') :
                IN_2 = gzip.open(fn2, 'rb')
            else :
                IN_2 = open(fn2, 'r')
        else :
            IN_2 = sys.stdin
    except IOError:
        sys.stderr.write("[Error] File %s cannot be open." % fn2)
        exit(-1)
    #
    line_1 = IN_1.readline()
    line_2 = IN_2.readline()
    #
    while line_1 and line_2 :
        try :
            chr_1, nuc_1, pos_1, pattern_1, dinuc_1, methyl_1, NmC_1, NC_1 = line_1.strip().split()
        except ValueError :
            sys.stderr.write("[Error] File [ %s ] may have wrong number of columns." % fn1)
            #exit(-1)
            line_1 = IN_1.readline()
            continue
        #
        try :
            chr_2, nuc_2, pos_2, pattern_2, dinuc_2, methyl_2, NmC_2, NC_2 = line_2.strip().split()
        except ValueError :
            sys.stderr.write("[Error] File [ %s ] may have wrong number of columns." % fn2)
            line_2 = IN_2.readline()
            continue
        #
        if chr_1 < chr_2 :
            line_1 = IN_1.readline()
        elif chr_1 > chr_2 :
            line_2 = IN_2.readline()
        else : # chr_1 == chr_2
            if int(pos_1) < int(pos_2) :
                line_1 = IN_1.readline()
            elif int(pos_1) > int(pos_2) :
                line_2 = IN_2.readline()
            else :
                if nuc_1 != nuc_2 or pattern_1 != pattern_2 or dinuc_1 != dinuc_2 :
                    sys.stderr.write("[Warning] Inconsistent information:")
                    sys.stderr.write("%s | %s" % (line_1, line_2) )
                #
                if CheckCtx( pattern_1, dinuc_1, ctx) :
                    print ("\t".join([chr_1, nuc_1, pos_1, pattern_1, dinuc_1,
                                 "%.2f" % float(methyl_1), NmC_1, NC_1,
                                 "%.2f" % float(methyl_2), NmC_2, NC_2 ] ) )
                               #  "%.2f" % float(methyl_1), "%d" % int(NmC_1), "%d" % int(NC_1),
                               #  "%.2f" % float(methyl_2), "%d" % int(NmC_2), "%d" % int(NC_2)]) )
                #
                line_1 = IN_1.readline()
                line_2 = IN_2.readline()
            #
        #
    #
    IN_1.close()
    if IN_2 is not sys.stdin:
        IN_2.close()
    #
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools intersect [-1 <CGmap_1>] -2 <CGmap_2> [-o <output>]\n" \
            "      (aka CGmapIntersect)\n" \
            "Description: \n" \
            "    Get the intersection of two CGmap files." \
            "Contact: Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-04-10\n" \
            "Output Format:\n" \
            "    Chr1  C  3541  CG  CG  0.8  4  5  0.4  4  10\n" \
            "When 1st CGmap file is:\n" \
            "    Chr1  C  3541  CG  CG  0.8  4  5\n" \
            ",and 2nd CGmap file is:\n" \
            "    Chr1  C  3541  CG  CG  0.4  4  10"
    #
    parser = OptionParser(usage)
    parser.add_option("-1", dest="CGmap_1", help="File name, end with .CGmap or .CGmap.gz. ", metavar="CGmap File")
    parser.add_option("-2", dest="CGmap_2", help="standard input if not specified", metavar="CGmap File")
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified. Compressed output if end with .gz")
    parser.add_option("-C", "--context", dest="context", default="",
                      help="specific context: CG, CH, CHG, CHH, CA, CC, CT, CW \
                            use all sites if not specified")
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
    CGmapIntersect(options.CGmap_1, options.CGmap_2, options.context)
#

# ===========================================
if __name__ == "__main__":
    main()
#
