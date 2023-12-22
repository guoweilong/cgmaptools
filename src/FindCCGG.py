#!/usr/bin/env python

"""
    cgmaptools - FindCCGG.py

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

#import os
import gzip
import sys
import re

def error(msg):
    sys.stderr.write( f'ERROR: {msg}')
    exit(1)
#

def read_fasta (fasta_file):
    """
        Iterates over all sequences in a fasta file. One at a time,
        without reading the whole file into the main memory.
    """
    #
    try:
        if fasta_file.endswith(".gz"):
            INPUT = gzip.open(fasta_file, "rt", encoding='UTF-8')
        else:
            INPUT = open(fasta_file, 'r')
    except IOError:
        print (f'[Error] Cannot find Fasta file : {fasta_file} !')
        exit(-1)
    sanitize = re.compile(r'[^ACTGN]')
    sanitize_seq_id = re.compile(r'[^A-Za-z0-9]')
    #
    chrome_seq = []
    chrome_id = None
    seen_ids = set()
    #
    for line in INPUT:
        if line[0] == '>':
            if chrome_id is not None:
                yield chrome_id, ''.join(chrome_seq)
            #
            chrome_id = sanitize_seq_id.sub('_', line.split()[0][1:])
            if chrome_id in seen_ids:
                error('Found identical sequence ids (id: %s) in the fasta file: %s.'
                      ' Please, make sure that all sequence ids are unique and '
                      ' contain only alphanumeric characters: A-Za-z0-9_' % (chrome_id, fasta_file))
            #
            seen_ids.add(chrome_id)
            chrome_seq = []
        else:
            chrome_seq.append(sanitize.sub('N', line.strip().upper()))
        #
    #
    yield chrome_id, ''.join(chrome_seq)
    INPUT.close()
#

def FindCCGG ( FastaFn, outputFn ):
    # the input fasta filename appended
    #
    cut_format = "CCGG"
    #
    total_chr = 0
    len_chr = dict()
    if outputFn:
        if outputFn.endswith("gz"):
            OUT = gzip.open(outputFn, 'wt', encoding='UTF-8')
        else:
            OUT = open(outputFn, 'w')
        #
    else:
        OUT = sys.stdout
    #
    for chr, seq in read_fasta(FastaFn):
        total_chr += 1
        L = len(seq)
        len_chr[chr] = L
        i = 1
        pre_pos = 0
        cur_pos = 0
        while i <= L - 4:
            if seq[i : i + 4] == "CCGG":
                cur_pos = i
                OUT.write("%s\n" % "\t".join( [chr, str(pre_pos+1), str(cur_pos+3)] ))
                pre_pos = cur_pos
            i += 1
        #
    #
#

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: cgmaptools findCCGG -i <genome.fa> [-o <output>]\n" \
            "      (aka FiindCCGG)\n" \
            "Description: Get the positions of all the C'CGG---CCG'G fragments.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last Update: 2018-01-02\n" \
            "Output Ex:\n" \
            "    chr1    4025    5652\n" \
            "    chr1    8274    8431"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Genome sequence file in Fasta format", metavar="FILE")
    parser.add_option("-o", dest="outfile",
                  help="Name of the output file (standard output if not specified). "
                       "Format: chr cCgg_pos ccGg_pos (0-base)\n", metavar="FILE")
    (options, args) = parser.parse_args()
    #
    if (options.infile is None):
        print(parser.print_help())
        exit(-1)
    #
    FindCCGG (options.infile, options.outfile)
#

# ===========================================
if __name__ == "__main__":
    main()
#
