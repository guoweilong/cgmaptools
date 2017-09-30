#!/usr/bin/env python

"""
    cgmaptools - ATCGmapStatCov.py

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

import gzip
import numpy as np

import datetime

from math import fsum
def average(x):
    return fsum(x)/float(len(x)) if x else 0
#
def logm(message):
    log_message = "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print log_message
#
def ATCGmapStatCov (fn, filetype = 'png', prefix = '', fW=8.0, fH=1.0):
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
        print "\n[Error]:\n\t File cannot be open: ", fn
        exit(-1)
    #
    line = IN.readline()
    CovDict = {}
    ChrLst = []
    chr_idx = -1
    ChrCovDict = []
    while line :
        try :
            chr, nuc, pos, pattern, dinuc, AW, TW, CW, GW, NW, AC, TC, CC, GC, NC, _ = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
            exit(-1)
        #
        COV = sum([int(i) for i in [AW, TW, CW, GW, NW, AC, TC, CC, GC, NC]])
        if COV in CovDict :
            CovDict[COV] += 1
        else :
            CovDict[COV] = 1
        #
        if chr not in ChrLst :
            ChrLst.append(chr)
            ChrCovDict.append({})
        #
        if COV in ChrCovDict[chr_idx]:
            ChrCovDict[chr_idx][COV] += 1
        else :
            ChrCovDict[chr_idx][COV] = 1
        #
        line = IN.readline()
        #
    #
    if IN is not sys.stdin:
        IN.close()
    #
    GlobalCov = (float(sum([i*j for i,j in CovDict.items()]))/sum([i for i in CovDict.values()]))
    print("OverAllCov\tglobal\t%.4f" % GlobalCov  )
    ChrCov = [0] * len(ChrLst)
    for chr_idx in xrange(len(ChrLst)):
        ChrCov[chr_idx] = float(sum([i * j for i, j in ChrCovDict[chr_idx].items()])) / sum([i for i in ChrCovDict[chr_idx].values()])
        print("OverAllCov\t%s\t%.4f" % (ChrLst[chr_idx], ChrCov[chr_idx] )   )
    for cov in sorted(CovDict.keys()) :
        print("CovAndCount\t%d\t%d" %(cov, CovDict[cov]) )
    #
    X = []
    Y = []
    for cov in sorted(CovDict.keys()):
        X.append(cov)
        Y.append(CovDict[cov])
        print("CovAndCount\t%d\t%d" % (cov, CovDict[cov]))
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
        plt.subplot(1, 2, 1)
        y_pos = range(len(ChrLst) + 1)
        plt.barh(y_pos, [GlobalCov] + ChrCov, align='center', alpha=0.6,
                 color=plt.cm.rainbow(np.linspace(1, 0, len(ChrLst) + 1)))
        plt.yticks(y_pos, ["Global"] + ChrLst)
        plt.title('Methylation Overall Coverage', fontsize=10)
        plt.xlim(0, max([GlobalCov] + ChrCov) * 1.2)
        # ==
        ax = plt.subplot(1, 2, 2)
        ax.plot(X, Y, '-o')
        ax.set_xscale("log")
        ax.yaxis.tick_right()
        plt.title('Distribution of Overall Coverage', fontsize=10)
        if prefix != "" :
            prefix = prefix + "."
        plt.savefig(prefix + "OverAllCove." + filetype, format=filetype)
        plt.clf()
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools oac stat [-i <ATCGmap>]\n" \
            "       (aka ATCGmapStatCov)\n" \
            "Description: Get the distribution of overall coverages.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com;\n" \
            "Last Update: 2016-12-16\n" \
            "Output Ex:\n" \
            "   OverAllCov      global  47.0395\n" \
            "   OverAllCov      chr1    45.3157\n" \
            "   OverAllCov      chr10   47.7380\n" \
            "   CovAndCount     1       1567\n" \
            "   CovAndCount     2       655\n" \
            "   CovAndCount     3       380"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="ATCGmap",
                      help= "File name end with .ATCGmap or .ATCGmap.gz. "
                            "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-f", "--figure-type", dest="FigType",
                      help="png, pdf, eps. Will not generate figure if not specified",
                      default=None, metavar="FILE")
    parser.add_option("-H", dest="fig_height", default=4, metavar="FLOAT",
                      help="Scale ratio for the Height of figure [Default: %default]")
    parser.add_option("-W", dest="fig_width",
                      help="Width of figure in inch [Default: %default]",
                      default=8, metavar="FLOAT")
    parser.add_option("-p", dest="prefix",
                      help="Prefix for output figures",
                      default="", metavar="STRING")  #
    #
    (options, args) = parser.parse_args()
    #
    #logm("Start analysis")
    ATCGmapStatCov(options.ATCGmap, options.FigType, options.prefix,
                   float(options.fig_width), float(options.fig_height))
    #logm("End of analysis")
#
# ===========================================
if __name__ == "__main__":
    main()
#
