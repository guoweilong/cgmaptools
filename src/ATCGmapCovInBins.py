#!/usr/bin/env python

"""
    cgmaptools - ATCGmapCovInBins.py

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


import gzip

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
def ATCGmapCovInBins (fn, step, filetype = 'png', prefix = '', title = '',
                       fH=4.0, fW=8.0):
    try:
        if fn :
            if fn.endswith(".gz") :
                IN = gzip.open(fn, 'rt',encoding='UTF-8')
            else :
                IN = open(fn, 'r')
            #
        else :
            IN = sys.stdin
        #
    except IOError:
        print ("\n[Error]:\n\t File cannot be open: %s" % fn)
        exit(-1)
    #
    line = IN.readline()
    posL = 1
    posR = step
    preChr = ""
    bin_list = []
    cov_list = []
    chr_list = []
    chr_start_list = []
    chr_end_list = []
    #
    while line :
        try :
            chr, nuc, pos, pattern, dinuc, AW, TW, CW, GW, NW, AC, TC, CC, GC, NC, _ = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
            exit(-1)
        #
        pos = int(pos)
        COV = sum([int(i) for i in [AW, TW, CW, GW, NW, AC, TC, CC, GC, NC]])
        if chr != preChr :
            if bin_list == [] :
                if preChr != "" :
                    print("%s\t%d\t%d\tna" % (preChr, posL, posR))
                    cov_list.append(float('nan'))
                #
            else :
                chr_end_list.append(len(cov_list))
                mean_cov = average(bin_list)
                print("%s\t%d\t%d\t%.4f" % (preChr, posL, posR, mean_cov ))
                cov_list.append(mean_cov)
                bin_list = []
            #
            chr_list.append( chr )
            chr_start_list.append( len(cov_list) )
            posL = 1
            posR = step
            preChr = chr
        #
        if pos <= posR :
            bin_list.append(COV)
        else :
            if bin_list == [] :
                print("%s\t%d\t%d\tna" % (preChr, posL, posR))
                cov_list.append(float('nan'))
            else :
                mean_cov = average(bin_list)
                print("%s\t%d\t%d\t%.4f" % (preChr, posL, posR, mean_cov))
                cov_list.append(mean_cov)
            #
            bin_list = []
            posL = posR+1
            posR += step
        #
        line = IN.readline()
    #
    if IN is not sys.stdin:
        IN.close()
    #
    if bin_list == []:
        if preChr != "":
            print("%s\t%d\t%d\tna" % (preChr, posL, posR))
            cov_list.append(float('nan'))
    else:
        mean_cov = average(bin_list)
        print("%s\t%d\t%d\t%.4f" % (preChr, posL, posR, mean_cov))
        cov_list.append(mean_cov)
    #
    chr_end_list.append( len(cov_list) )
    #
    if filetype in ['png', 'eps', 'pdf'] :
        import matplotlib
        # Force matplotlib to not use any Xwindows backend.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcdefaults()
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import cm
        # ===
        chr_mid_list  = [ (s+e)/2 for s,e in zip(chr_start_list, chr_end_list)]
        plt.figure(figsize=(fW, fH))
        MaxCOV = NanMax(cov_list) * 1.1 + NanMin(cov_list)
        plt.ylim([0, MaxCOV])
        plt.plot(range(0, len(cov_list), 1), cov_list, 'b-')
        for x in chr_start_list[1:] :
            plt.plot( [x, x], [-MaxCOV, MaxCOV], 'k-' )
        #
        plt.xticks(chr_mid_list, chr_list, rotation=25)
        plt.ylabel("Overall coverage")
        plt.title(title)
        #
        plt.savefig(prefix + ".OverallCovInBins."+filetype, format=filetype)
        plt.clf()
    #
#
#
from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools oac bin  [-i <ATCGmap>] [-B 5000000] \n" \
            "       (aka ATCGmapCovInBins)\n"\
            "Description: Generate the overall coverage in Bins.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com;\n" \
            "Last Update: 2016-12-07\n" \
            "Output Ex:\n" \
            "   chr1    1       5000    29.0000\n" \
            "   chr1    5001    10000   30.0396\n" \
            "   chr2    1       5000    35.0755\n" \
            "   chr2    5001    10000   40.0027\n" \
            "   chr3    1       5000    na"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="ATCGmap",
                      help="File name end with .ATCGmap or .ATCGmap.gz. "
                           "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-B", dest="bin_size", default=5000000,
                      help="Define the size of bins [Default: %default]")
    parser.add_option("-f", "--figure-type", dest="FigType",
                      help="png, pdf, eps. Will not generate figure if not specified",
                      default=None, metavar="FILE")
    parser.add_option("-H", dest="fig_height", default=4, metavar="FLOAT",
                      help="Height of figure in inch [Default: %default]")
    parser.add_option("-W", dest="fig_width",
                      help="Width of figure in inch [Default: %default]",
                      default=8, metavar="FLOAT")
    parser.add_option("-p", dest="prefix",
                      help="Prefix for output figures",
                      default="", metavar="STRING")
    parser.add_option("-t", "--title", dest="title",
                      help="title in the output figures",
                      default="", metavar="STRING")
    #
    (options, args) = parser.parse_args()
    #
    ATCGmapCovInBins(options.ATCGmap, int(options.bin_size), options.FigType, options.prefix,
                     options.title, float(options.fig_height), float(options.fig_width))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
