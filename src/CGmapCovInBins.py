#!/usr/bin/env python

"""
    cgmaptools - CGmapCovInBins.py

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
def CGmapCovInBins (fn, step, CTX, filetype = 'png', prefix = '', title = '', fH = 4.0, fW = 8.0):
    try:
        if fn:
            if fn.endswith(".gz"):
                IN = gzip.open(fn, 'rt', encoding='UTF-8')
            else:
                IN = open(fn, 'r')
            #
        else:
            IN = sys.stdin
        #
    except IOError:
        print(f'\n[Error]:\n\t File cannot be open: {fn}')
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
    while line:
        try:
            chr, nuc, pos, pattern, dinuc, methyl, MC, NC = line.strip().split()
        except ValueError:
            print(f'\n[Error]:\n\t File [ {fn} ] may have wrong number of columns.')
            exit(-1)
        #
        pos = int(pos)
        NC = float(NC)
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
        if chr != preChr:
            if bin_list == [] :
                if preChr != "":
                    print(f'{preChr}\t{posL}\t{posR}\tna')
                    cov_list.append(float('nan'))
                #
            else:
                chr_end_list.append(len(cov_list))
                mean_cov = average(bin_list)
                print(f'{preChr}\t{posL}\t{posR}\t{mean_cov:.4f}')
                cov_list.append(mean_cov)
                bin_list = []
            #
            chr_list.append(chr)
            chr_start_list.append(len(cov_list))
            posL = 1
            posR = step
            preChr = chr
        #
        if pos <= posR:
            bin_list.append(float(NC))
        else:
            if bin_list == []:
                print(f'{preChr}\t{posL}\t{posR}\tna')
                cov_list.append(float('nan'))
            else:
                mean_cov = average(bin_list)
                print(f'{preChr}\t{posL}\t{posR}\t{mean_cov:.4f}')
                cov_list.append(mean_cov)
            #
            bin_list = []
            posL = posR+1
            posR += step
        #
        line = IN.readline()
        #
    #
    if bin_list == []:
        if preChr != "":
            print(f'{preChr}\t{posL}\t{posR}\tna')
            cov_list.append(float('nan'))
        #
    else:
        mean_cov = average(bin_list)
        print(f'{preChr}\t{posL}\t{posR}\t{mean_cov:.4f}')
        cov_list.append(mean_cov)
    #
    chr_end_list.append(len(cov_list))
    if IN is not sys.stdin:
        IN.close()
    #
    if filetype in ['png', 'eps', 'pdf']:
        import matplotlib
        # Force matplotlib to not use any Xwindows backend.
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcdefaults()
        import matplotlib.pyplot as plt
        # ===
        chr_mid_list  = [(s+e)/2 for s,e in zip(chr_start_list, chr_end_list)]
        plt.figure(figsize=(fW, fH))
        MaxCOV = NanMax(cov_list) * 1.1 + NanMin(cov_list)
        plt.ylim([0, MaxCOV])
        plt.plot(range(0, len(cov_list), 1), cov_list, 'b-')
        for x in chr_start_list[1:] :
            plt.plot( [x, x], [-MaxCOV, MaxCOV], 'k-' )
        #
        plt.xticks(chr_mid_list, chr_list, rotation=25)
        plt.ylabel("Methylaiton effective coverage")
        plt.title(title)
        #
        if prefix != "":
            prefix = prefix + "."
        plt.savefig(prefix + "MethEffeCovInBins."+filetype, format=filetype)
        plt.clf()
    #
#
#
from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools mec bin [-i <CGmap>] [-B 5000000] \n" \
            "      (aka CGmapCovInBins)\n" \
            "Description: Generate the methylation-effective coverage in Bins.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com;\n" \
            "Last Update: 2018-01-02\n" \
            "Output Ex:\n" \
            "   chr1    1       5000    29.0000\n" \
            "   chr1    5001    10000   30.0396\n" \
            "   chr2    1       5000    35.0755\n" \
            "   chr2    5001    10000   40.0027\n" \
            "   chr3    1       5000    na"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap",
                      help= "File name end with .CGmap or .CGmap.gz. "
                            "If not specified, STDIN will be used.", metavar="FILE")
    #
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
    parser.add_option("-C", "--context", dest="context", default="",
                      help="specific context: CG, CH, CHG, CHH, CA, CC, CT, CW \
                            use all sites if not specified")
    #
    (options, args) = parser.parse_args()
    #
    CGmapCovInBins(options.CGmap, int(options.bin_size), options.context, options.FigType, options.prefix, options.title,
                   float(options.fig_height), float(options.fig_width))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
