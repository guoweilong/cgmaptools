#!/usr/bin/env python

"""
    cgmaptools - CGmapsMethInBins.py

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

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

import gzip


from math import fsum
from math import isnan

def average(x):
    return fsum(x)/float(len(x)) if x else 0
#
def NanMax(x):
    return max([i for i in x if not isnan(i)])
#
def NanMin(x):
    return min([i for i in x if not isnan(i)])
#

def CGmapMethylInBins (fn_lst_str, tag_lst_str, coverage, coverageXY, step, CTX = ""):
    fn_lst = fn_lst_str.split(",")
    N_file = len(fn_lst)
    chr_lst = []
    chr_idx = {}
    meth_DB = []
    tag_lst = tag_lst_str.split(",")
    if (len(tag_lst)!= len(fn_lst)) :
        print("[Error]: tag list and file list are not same in lengths.")
        exit(-1)
    #
    #
    # works for python2 and python3
    try:
        xrange
    except NameError:
        xrange = range
    #
    for fn_id in xrange( len(fn_lst) ) :
        fn = fn_lst[fn_id]
        try:
            if fn :
                if fn.endswith(".gz") :
                    IN = gzip.open(fn, 'rb')
                else :
                    IN = open(fn, 'r')
                #
            else :
                print("[Error]:\n\t un-specified name")
                exit(-1)
            #
        except IOError:
            print("\n[Error]:\n\t File cannot be open: %s" % fn)
            exit(-1)
        #
        line = IN.readline()
        posL = 1
        posR = step
        preChr = ""
        bin_list = []
        #
        reg_id = 0
        while line :
            try :
                chr, nuc, pos, pattern, dinuc, methyl, MC, NC = line.strip().split()
            except ValueError :
                print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
                exit(-1)
            #
            if CTX not in ["", "C"]:
                if CTX in ["CG", "CHG", "CHH"]:
                    if pattern != CTX:
                        line = IN.readline()
                        continue
                        #
                elif CTX in ["CA", "CC", "CT"]:
                    if dinuc != CTX:
                        line = IN.readline()
                        continue
                        #
                elif CTX == "CH":
                    if dinuc not in ["CHG", "CHH"]:
                        line = IN.readline()
                        continue
                        #
                elif CTX == "CW":
                    if dinuc not in ["CA", "CC"]:
                        line = IN.readline()
                        continue
                        #
                else:
                    line = IN.readline()
                    continue
                    #
            #
            pos = int(pos)
            MC = int(MC)
            NC = int(NC)
            methyl = float(MC)/NC
            if chr not in chr_lst:
                chr_idx[chr] = len(chr_lst)
                chr_lst.append(chr)
                meth_DB.append([])
            #
            if  (NC>=coverage) or ((chr in ["chrX", "chrY", "ChrX", "ChrY"]) and (NC>=coverageXY)):
                if chr != preChr : # new chr
                    if preChr != "" :
                        if reg_id >= len(meth_DB[ chr_idx[preChr] ]) :
                            meth_DB[chr_idx[preChr]].append([float('nan') for i in range(N_file)])
                        #
                        if bin_list != [] :
                            meth_DB[chr_idx[preChr]][reg_id][fn_id] = average(bin_list)
                            bin_list = []
                        #
                    #
                    posL = 1
                    posR = step
                    preChr = chr
                    reg_id = 0
                #
                while posR < pos :
                    if reg_id >= len(meth_DB[ chr_idx[preChr] ]) :
                        meth_DB[chr_idx[preChr]].append([float('nan') for i in range(N_file)])
                    #
                    meth_DB[chr_idx[preChr]][reg_id][fn_id] = average(bin_list)
                    bin_list = []
                    posL = posR+1
                    posR += step
                    reg_id += 1
                #
                bin_list.append(methyl)
                #
            #
            line = IN.readline()
        #
        if preChr != "":
            if reg_id >= len(meth_DB[chr_idx[preChr]]):
                # print "===5==="
                meth_DB[chr_idx[preChr]].append([float('nan') for i in range(N_file)])
            #
            if bin_list != []:
                meth_DB[chr_idx[preChr]][reg_id][fn_id] = average(bin_list)
                bin_list = []
            #
        #
        if IN is not sys.stdin:
            IN.close()
        #
    #
    print("\t".join(["chr", "start", "end"] + tag_lst))
    for chr_idx in xrange(len(chr_lst)) :
        chr = chr_lst[chr_idx]
        for reg_idx in xrange(len(meth_DB[chr_idx])) :
            #print meth_DB
            meth_vec = meth_DB[chr_idx][reg_idx]
            #print meth_vec
            print("%s\t%d\t%d\t%s" % (chr, reg_idx*step+1, (reg_idx+1)*step,
                                      '\t'.join([ ("%.2f"%meth) for meth in meth_vec]) ))
            #
        #
    #
#
from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools mmbin [-l <1.CGmap[,2.CGmap,..]>] [-c 10 --CXY 5 -B 5000000] \n" \
            "      (aka CGmapsMethInBins)\n" \
            "Description: Generate the methylation in Bins.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-05-02\n" \
            "Output Ex:\n" \
            "   chr1    1       5000    0.0000\n" \
            "   chr1    5001    10000   0.0396\n" \
            "   chr2    1       5000    0.0755\n" \
            "   chr2    5001    10000   0.0027\n" \
            "   chr3    1       5000    na"
    #
    parser = OptionParser(usage)
    parser.add_option("-l", dest="CGmapLst", help="File name list, end with .CGmap or .CGmap.gz. "
                                               "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-t", dest="TagLst", help="List of samples", default = None, metavar="FILE")

    parser.add_option("-B", dest="bin_size", default=5000000,
                      help="Define the size of bins [Default: %default]")
    parser.add_option("-C", "--context", dest="context", default="",
                      help="specific context: CG, CH, CHG, CHH, CA, CC, CT, CW \
                            use all sites if not specified")
    parser.add_option("-c", dest="coverage", default=10,
                      help="The minimum coverage for site selection [Default: %default]")
    parser.add_option("--cXY", dest="coverageXY", default=None,
                      help="Coverage for chrX/Y should be half that of autosome for male [Default: same with -c]")
    #
    (options, args) = parser.parse_args()
    if options.coverageXY is None :
        options.coverageXY = options.coverage
    if options.TagLst is None :
        options.TagLst = options.CGmapLst
    #
    CGmapMethylInBins(options.CGmapLst, options.TagLst, int(options.coverage), int(options.coverageXY),
                      int(options.bin_size), options.context )
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
