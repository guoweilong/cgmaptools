#!/usr/bin/env python

"""
    cgmaptools - CGmapStatMeth.py

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

"""
Chr1    C       3541    CG      CG      0.0     0       1
Chr1    C       3548    CHH     CC      0.0     0       1
Chr1    C       3549    CHH     CA      0.0     0       1
"""

import gzip
import numpy as np

import datetime

def logm(message):
    log_message = "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print (log_message)
#

def GetQuant( methyl, NQuant = 5 ) :
    return int(np.ceil(methyl * NQuant))
#

def CGmapStatMeth (fn, coverage = 10, filetype = "png", prefix = "", title="", fH=4.0, fW=8.0):
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
        print ("\n[Error]:\n\t File cannot be open: %s" % fn)
        exit(-1)
    #
    #print "Read in the file"
    line = IN.readline()
    #               0     1     2     3      4     5      6    7     8
    context_lst = ['C', 'CG', 'CHG', 'CHH', 'CA', 'CC', 'CT', 'CH', 'CW']
    context_idx = {'C':0, 'CG':1, 'CHG':2, 'CHH':3, 'CA':4, 'CC':5, 'CT':6, 'CH':7, 'CW':8}
    chr_lst = []
    sum_mC_byChr   = {}
    count_mC_byChr = {}
    NQuant = 5
    quant_mC = [ [0]*(NQuant+1) for i in xrange(len(context_lst)) ]
    #
    # ==================
    # Read in file by line
    while line :
        try :
            chr, nuc, pos, pattern, dinuc, methyl, MC, NC = line.strip().split()
        except ValueError :
            print("\n[Error]:\n\t File [ %s ] may have wrong number of columns." % fn)
            exit(-1)
        #
        MC = int(MC)
        NC = int(NC)
        if NC < coverage or NC <= 0 :
            line = IN.readline()
            continue
        methyl = float(MC)/NC
        if chr not in chr_lst :
            chr_lst.append(chr)
            sum_mC_byChr[chr]   = [0.0] * len(context_lst)
            count_mC_byChr[chr] = [0] * len(context_lst)
        #
        qt = GetQuant(methyl)
        #
        sum_mC_byChr[chr][context_idx['C']] += methyl
        count_mC_byChr[chr][context_idx['C']] += 1
        quant_mC[context_idx['C']][qt] += 1
        if pattern == "CG" :
            sum_mC_byChr[chr][context_idx['CG']] += methyl
            count_mC_byChr[chr][context_idx['CG']] += 1
            quant_mC[context_idx['CG']][qt] += 1
        else:
            sum_mC_byChr[chr][context_idx['CH']] += methyl
            count_mC_byChr[chr][context_idx['CH']] += 1
            quant_mC[context_idx['CH']][qt] += 1
            if pattern == "CHG" :
                sum_mC_byChr[chr][context_idx['CHG']] += methyl
                count_mC_byChr[chr][context_idx['CHG']] += 1
                quant_mC[context_idx['CHG']][qt] += 1
            elif pattern == "CHH" :
                sum_mC_byChr[chr][context_idx['CHH']] += methyl
                count_mC_byChr[chr][context_idx['CHH']] += 1
                quant_mC[context_idx['CHH']][qt] += 1
            #
            if dinuc == "CA" :
                sum_mC_byChr[chr][context_idx['CA']] += methyl
                count_mC_byChr[chr][context_idx['CA']] += 1
                quant_mC[context_idx['CA']][qt] += 1
                #
                sum_mC_byChr[chr][context_idx['CW']] += methyl
                count_mC_byChr[chr][context_idx['CW']] += 1
                quant_mC[context_idx['CW']][qt] += 1
            elif dinuc == "CC" :
                sum_mC_byChr[chr][context_idx['CC']] += methyl
                count_mC_byChr[chr][context_idx['CC']] += 1
                quant_mC[context_idx['CC']][qt] += 1
            elif dinuc == "CT" :
                sum_mC_byChr[chr][context_idx['CT']] += methyl
                count_mC_byChr[chr][context_idx['CT']] += 1
                quant_mC[context_idx['CT']][qt] += 1
                #
                sum_mC_byChr[chr][context_idx['CW']] += methyl
                count_mC_byChr[chr][context_idx['CW']] += 1
                quant_mC[context_idx['CW']][qt] += 1
            #
        #
        #print quant_mC
        line = IN.readline()
    #
    #print quant_mC
    if IN is not sys.stdin:
        IN.close()
    # =======================
    # Summary Section
    mean_mC = [0.0] * len(context_lst)
    count_C = [0] * len(context_lst)
    sd_mC = [np.NaN] * len(context_lst)
    for i in xrange( len(context_lst) ) :
        Denominator = sum([count_mC_byChr[chr][i] for chr in count_mC_byChr])
        if Denominator not in [0, float('nan')] :
            mean_mC[i] = sum([sum_mC_byChr[chr][i] for chr in sum_mC_byChr])/Denominator
        else :
            mean_mC[i] = float('nan')
        #
        if len(chr_lst) > 2 :
            sd_mC[i] = np.std( [ sum_mC_byChr[chr][i]/count_mC_byChr[chr][i] if count_mC_byChr[chr][i]>0 else np.NaN
                                 for chr in sum_mC_byChr ] )
        count_C[i] = sum( count_mC_byChr[chr][i] for chr in chr_lst)
    #
    contrib_mC = [i*j for i,j in zip(mean_mC, count_C)]
    contrib_mC_sum = contrib_mC[0]
    if contrib_mC_sum not in [0, float('nan')] :
        contrib_mC = [i/contrib_mC_sum for i in contrib_mC]
    else :
        contrib_mC = [float('nan') for i in contrib_mC]
    #
    # -----------------
    # stat summary for BulkMethyl
    print ( "\t".join(["MethStat", "context"] + context_lst) )
    print ( "\t".join(["mean_mC ", "global"] + ["%.4f"%i for i in mean_mC]) )
    print ( "\t".join(["sd_mCbyChr", "global"] + ["%.4f"%i for i in sd_mC]) )
    print ( "\t".join(["count_C ", "global"] + ["%d"%i for i in count_C]) )
    print ( "\t".join(["contrib_mC ", "global"] + ["%.4f"%i for i in contrib_mC]) )
    print ( "\t".join(["quant_mC ", "[0]"] + ["%d" % j[0] for j in quant_mC]) )
    for i in xrange(NQuant) :
        print ("\t".join(["quant_mC ", "(%.2f,%.2f]"%(float(i)/NQuant, float(i+1)/NQuant)] + ["%d"%(j[i+1]) for j in quant_mC]) )
        #
    for chr in chr_lst :
        print ("\t".join(["mean_mC_byChr", chr] + ["%.4f"%(sum_mC_byChr[chr][i]/count_mC_byChr[chr][i])
                                        if count_mC_byChr[chr][i]>0 else "NaN" for i in xrange( len(context_lst) ) ]) )
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
        plt.figure(figsize=(fW, fH))
        x_pos = np.arange(len(context_lst))
        plt.bar(x_pos, mean_mC, yerr = sd_mC, align = 'center', alpha = 0.6,
                color = plt.cm.rainbow(np.linspace(1,0, len(context_lst) ) ) )
        plt.xticks(x_pos, context_lst)
        plt.ylabel('Average methylation level')
        plt.ylim([0,1])
        plt.title(title)
        #plt.legend()
        if prefix != "" :
            prefix = prefix + "."
        plt.savefig(prefix + "BulkMeth."+filetype, format=filetype)
        plt.clf()
        # ===
        plt.figure(figsize=(fW, fH))
        NContext = len(context_lst)
        for i in xrange( NContext ) :
            plt.subplot(1, NContext, i+1)
            y_pos = np.arange(NQuant+1)
            plt.barh(y_pos, quant_mC[i], align = 'center', alpha = 0.6,
                color = plt.cm.rainbow(1.0-float(i)/NContext), linewidth = 0 )
            plt.autoscale(enable=True, axis='x')
            plt.xticks([])
            if i == 0 :
                plt.yticks(y_pos, ['[0]'] + ["(%.2f, %.2f]"%(float(j)/NQuant, float(j+1)/NQuant) for j in xrange(NQuant)])
            else :
                plt.yticks(y_pos,[])
            #
            plt.title(context_lst[i])
        #
        plt.savefig(prefix + "MethFrag."+filetype, format=filetype)
        plt.clf()
        # ===
        plt.figure(figsize=(fW, fH*0.8))
        j = 1
        for CONTEXT in [ ["CG", "CHG", "CHH"], ["CG", "CW", "CC"], ["CG", "CA", "CC", "CT"]] :
            plt.subplot(1, 3, j)
            Freq = [contrib_mC[context_idx[i]] for i in CONTEXT]
            FreqSum = sum(Freq)
            Freq = [i/FreqSum for i in Freq]
            #print Freq
            ColPanel = [plt.cm.rainbow(1.0 - float(context_idx[i]) / NContext, alpha = 0.6) for i in CONTEXT]
            plt.pie(Freq, labels = CONTEXT, startangle=90, colors = ColPanel )
            j = j + 1
        #
        plt.savefig(prefix + "MethContri."+filetype, format=filetype)
        plt.clf()
    #
#
#
from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: cgmaptools mstat [-i <CGmap>]  \n" \
            "      (aka CGmapStatMeth)\n" \
            "Description: Generate the bulk methylation.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-01-02\n" \
            "Output Ex:\n" \
            "   MethStat        context C       CG      CHG	    CHH	    CA      CC      CT      CH      CW\n" \
            "   mean_mC         global  0.0798  0.3719  0.0465  0.0403  0.0891  0.0071  0.0241  0.0419  0.0559\n" \
            "   sd_mCbyChr      global  0.0078  0.0341  0.0163  0.0110  0.0252  0.0049  0.0076  0.0096  0.0148\n" \
            "   count_C         global  10000   1147    2332    6521    3090    2539    3224    8853    6314\n" \
            "   contrib_mC      global  1.0000  0.5348  0.1360  0.3292  0.3452  0.0228  0.0973  0.4652  0.4424\n" \
            "   quant_mC        [0]     8266    471     2012    5783    2422    2421    2952    7795    5374\n" \
            "   quant_mC   (0.00 ,0.20] 705     182     155     368     272     97      154     523     426\n" \
            "   mean_mC_byChr   chr1    0.0840  0.4181  0.0340  0.0412  0.0794  0.0065  0.0251  0.0393  0.0513\n" \
            "   mean_mC_byChr   chr10   0.0917  0.4106  0.0758  0.0421  0.0968  0.0097  0.0349  0.0502  0.0655"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap",
                      help="File name end with .CGmap or .CGmap.gz. "
                           "If not specified, STDIN will be used.", metavar="FILE")
    parser.add_option("-c", dest="coverage", default=10,
                      help="The minimum coverage for site selection [Default: %default]")
    parser.add_option("-f", "--figure-type", dest="FigType",
                      help="png, pdf, eps. Will not generate figure if not specified",
                      default = None, metavar="FILE")
    parser.add_option("-H", dest="fig_height", default=3, metavar="FLOAT",
                      help="Height of figure in inch [Default: %default]")
    parser.add_option("-W", dest="fig_width",
                      help="Width of figure in inch [Default: %default]",
                      default=8, metavar="FLOAT")
    parser.add_option("-p", dest="prefix",
                      help="Prefix for output figures",
                      default = "", metavar="STRING")
    parser.add_option("-t", "--title", dest="title",
                      help="title in the output figures",
                      default="", metavar="STRING")
    #
    (options, args) = parser.parse_args()
    CGmapStatMeth(options.CGmap, int(options.coverage), options.FigType,
                  options.prefix, options.title,
                   float(options.fig_height), float(options.fig_width))
    #
#

# ===========================================
if __name__ == "__main__":
    main()
#
