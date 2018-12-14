#!/usr/bin/env python

"""
    cgmaptools - CGmapSplitByChr.py

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
import numpy as np

import datetime

from math import fsum
from math import isnan
#
def average(x):
    return fsum(x)/float(len(x)) if x else 0
#
def NanMax(x):
    return max([i for i in x if not isnan(i)])
#
def NanMin(x):
    return min([i for i in x if not isnan(i)])
#
def revcumsum(seq):
    cumseq = []
    s = 0
    for c in reversed(seq):
       s+= c
       cumseq.append(s)
    return list(reversed(cumseq))
#
def logm(message):
    log_message = "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print(log_message)

def CGmapStatCov (fn, CTX, filetype = 'png', prefix = '', fH=1.0, fW=11.0):
    try:
        if fn :
            if fn.endswith(".gz") :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
        else :
            IN = sys.stdin
    except IOError:
        print ("\n[Error]:\n\t File cannot be open: %s" % fn)
        exit(-1)
    #
    line = IN.readline()
    CovDict = {}
    ChrLst = []
    chr_idx = -1
    ChrCovDict = []
    while line :
        try :
            chr, nuc, pos, pattern, dinuc, methyl, MC, NC = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
            exit(-1)
        #
        # check if line could be filtered by context
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
        NC = int(NC)
        if NC in CovDict :
            CovDict[NC] += 1
        else :
            CovDict[NC] = 1
        #
        if chr not in ChrLst :
            ChrLst.append(chr)
            ChrCovDict.append({})
        #
        if NC in ChrCovDict[chr_idx]:
            ChrCovDict[chr_idx][NC] += 1
        else :
            ChrCovDict[chr_idx][NC] = 1
        #
        line = IN.readline()
    #
    if IN is not sys.stdin:
        IN.close()
    #
    GlobalCovSum = sum([i for i in CovDict.values()])
    if GlobalCovSum > 0 :
        GlobalCov = (float(sum([i*j for i,j in CovDict.items()]))/GlobalCovSum)
        print("MethEffectCov\tglobal\t%.4f" % GlobalCov)
    else :
        GlobalCovSum = 0
        print("MethEffectCov\tglobal\tNaN")
    #
    #
    # works for python2 and python3
    try:
        xrange
    except NameError:
        xrange = range
    #
    ChrCov = [0] * len(ChrLst)
    for chr_idx in xrange(len(ChrLst)):
        ChrCov[chr_idx] = float(sum([i * j for i, j in ChrCovDict[chr_idx].items()])) / sum([i for i in ChrCovDict[chr_idx].values()])
        print("MethEffectCov\t%s\t%.4f" % (ChrLst[chr_idx], ChrCov[chr_idx]))
    X = []
    Y = []
    for cov in sorted(CovDict.keys()) :
        X.append(cov)
        Y.append(CovDict[cov])
        print("CovAndCount\t%d\t%d" %(cov, CovDict[cov]) )
    #
    if filetype in ['png', 'eps', 'pdf']:
        import matplotlib
        # Force matplotlib to not use any Xwindows backend.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcdefaults()
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import cm
        # ===
        plt.figure(figsize=(fW, len(ChrLst)*0.22*fH))
        plt.subplot(1, 3, 1)
        y_pos = range(len(ChrLst)+1)
        plt.barh(y_pos, [GlobalCov] + ChrCov, align='center', alpha=0.6,
                color=plt.cm.rainbow(np.linspace(1, 0, len(ChrLst)+1)))
        plt.yticks(y_pos, ["Global"]+ChrLst)
        plt.title('Methylation Effective Coverage', fontsize=10)
        plt.xlim(0, max([GlobalCov] + ChrCov)*1.2)
        # ==
        ax = plt.subplot(1, 3, 2)
        plt.plot(X, Y, '-o')
        ax.yaxis.tick_right()
        ax.set_xscale("log")
        plt.title('Distribution of Effective Coverage', fontsize=10)
        # ==
        ax = plt.subplot(1, 3, 3)
        Z = revcumsum(Y)
        plt.plot(X, Z, '-o')
        ax.yaxis.tick_right()
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.title('Cumulative Effective Coverage', fontsize=10)
        if prefix != "" :
            prefix = prefix + "."
        #
        plt.savefig(prefix + "MethEffectCove." + filetype, format=filetype)
        plt.clf()
        # ===

        # ===
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools mec stat [-i <CGmap>]  \n" \
            "      (aka CGmapStatCov)\n" \
            "Description: Get the distribution of methylation-effective coverages.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-05-02\n" \
            "Output Ex:\n" \
            "   MethEffectCove  global  47.0395\n" \
            "   MethEffectCove  chr1    45.3157\n" \
            "   MethEffectCove  chr10   47.7380\n" \
            "   CovAndCount     1       1567\n" \
            "   CovAndCount     2       655\n" \
            "   CovAndCount     3       380"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap",
                      help="File name end with .CGmap or .CGmap.gz. "
                            "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-f", "--figure-type", dest="FigType",
                      help="png, pdf, eps. Will not generate figure if not specified",
                      default=None, metavar="FILE")
    parser.add_option("-H", dest="fig_height", default=4, metavar="FLOAT",
                      help="Scale factor for the Height of figure [Default: %default]")
    parser.add_option("-W", dest="fig_width",
                      help="Width of figure in inch [Default: %default]",
                      default=11, metavar="FLOAT")
    parser.add_option("-p", dest="prefix", help="Prefix for output figures",
                      default="", metavar="STRING")  #
    parser.add_option("-C", "--context", dest="context", default="",
                      help="specific context: CG, CH, CHG, CHH, CA, CC, CT, CW \
                            use all sites if not specified")
    (options, args) = parser.parse_args()
    #
    #logm("Start analysis")
    CGmapStatCov(options.CGmap, options.context, options.FigType, options.prefix,
                 float(options.fig_height), float(options.fig_width))
    #logm("End of analysis")
#
# ===========================================
if __name__ == "__main__":
    main()
#
