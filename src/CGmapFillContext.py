#!/usr/bin/env python

"""
    cgmaptools - CGmapFIllContext.py

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

#import math
#x=float('nan')
#math.isnan(x)

import gzip

def error(msg):
    sys.stderr.write( "ERROR: %s" % msg )
    exit(1)

genome=dict()
genome_len=dict()

def read_genome (fasta_file) :
    """
        Iterates over all sequences in a fasta file. One at a time,
        without reading the whole file into the main memory.
    """
    #
    try :
        INPUT = (gzip.open if fasta_file.endswith('.gz') else open)(fasta_file)
    except IOError:
        print ("[Error] Cannot find Fasta file : %s !" % fasta_file)
        exit(-1)
    #
    chr = ""
    seq_list = []
    for line in INPUT :
        if line[0] == '>':
            if chr != "" and seq_list != [] :
                genome[chr]="".join(seq_list)
                seq_list=[]
            #
            chr = line[1:].strip()
            sys.stderr.write( "# Reading " + chr )
        else:
            seq_list.append(line.strip().upper())
        #
    #
    if chr != "" and seq_list != [] :
        genome[chr] = "".join(seq_list)
        seq_list=[]
    #
    INPUT.close()
    sys.stderr.write("# Finished reading the genome")
    for chr in genome :
        genome_len[chr] = len(genome[chr])
        sys.stderr.write( "# %s\t%d" % (chr, genome_len[chr]) )
    #
#
bc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
           'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
#
def bc ( ch ) :
    if ch in bc_dict :
        return bc_dict[ch]
    else :
        return ch
    #
#

def AntisenseRead ( read ) :
    read="".join([ bc(i) for i in read])
    return read[::-1]
#

def CGmapFillContext(CGmapF, genomeF, base=1) :
    read_genome(genomeF)
    # Read all the CGmap files
    try :
        if CGmapF :
            if CGmapF.endswith('.gz') :
                IN = gzip.open(CGmapF, 'rb')
            else :
                IN = open(CGmapF, 'r')
            #
        else :
            IN = sys.stdin
        #
    except IOError :
        print ("\n[Error]:\n\t File cannot be open: %s " % CGmapF )
        exit(-1)
    #
    # chr1    C       4654    CG      CG      0.84  11      13
    tokens = []
    for line in IN :
        tokens = line.strip().split()
        chr = tokens[0]
        nuc = tokens[1]
        pos = int(tokens[2]) - base
        if 2< pos < genome_len[chr] :
            if nuc == "C" :
                trimer=genome[chr][pos:(pos+3)]
            elif nuc == "G" :
                trimer=AntisenseRead(genome[chr][(pos-2):(pos+1)])
            else :
                continue
                # warning : not CGmap file
            #
            #print trimer
            dimer=trimer[0:2]
            #print trimer + "\t" + dimer
            if trimer[0] != "C" :
                sys.stderr.write( "# [warning] Wrong NUC: %s\t%s\t%s\t%s\n" % (chr, nuc, pos, trimer) )
                continue
            #
            try :
                if dimer == "CG" :
                    pattern = "CG"
                elif trimer[2] == "G" :
                    pattern = "CHG"
                else :
                    pattern = "CHH"
                #
                tokens[3] = pattern
                tokens[4] = dimer
            except IndexError :
                continue
            #
            print ("\t".join(tokens) )
        else :
            print ("\t".join(tokens) )
        #
    #
    IN.close()
#


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools refill [-i <CGmap>] -g <genome.fa> [-o output]\n" \
            "      (aka CGmapFillContext)\n" \
            "Description: Fill the CG/CHG/CHH and CA/CC/CT/CG context.\n" \
            "             Other fields will not be affected.\n" \
            "             Can be applied to ATCGmap file.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com; \n" \
            "Last Update: 2018-01-02\n" \
            "Index Ex:\n" \
            "   Chr1    C       3541    -       -       0.0     0       1\n" \
            "Output Ex:\n" \
            "   Chr1    C       3541    CG      CG      0.0     0       1"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="CGmap",
                      help="Input CGmap file (CGmap or CGmap.gz)", metavar="STRING", default=None)
    parser.add_option("-g", dest="genome",
                      help="genome file, FASTA format (gzipped if end with \'.gz\')", metavar="STRING")
    parser.add_option("-o", dest="outfile",
                      help="Output file name (gzipped if end with \'.gz\')", metavar="STRING")
    parser.add_option("-0", "--0-base", dest="base0", action="store_true",
                      help="0-based genome if specified [Default: 1-based]", default=False)
    (options, args) = parser.parse_args()
    #
    if ( options.genome is None) :
        parser.print_help()
        exit(-1)
    #
    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')
        #
    if options.base0 :
        base = 0
    else :
        base = 1
    #
    CGmapFillContext(options.CGmap, options.genome, base)
#

# ===========================================
if __name__ == "__main__":
    main()
#
