#!/usr/bin/env python

"""
    cgmaptools - SNVFromATCGmap.py

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

# Guo, Weilong; guoweilong@126.com; created on 2014-12-15

import sys
#import os
#import os.path
#import re


def binomialCoeff(n, k):
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result

def dbinom(n, k, prob):
    return binomialCoeff(n, k) * (prob ** k) * ((1-prob) ** (n-k))

def qbinom(p, size, prob, lower_tail = True):
    '''Binomial quantiles'''
    if not lower_tail :
        prob = 1 - prob
    sum = 0
    for i in xrange(size) :
        p_binom = dbinom(size, i, prob)
        if (sum + p_binom > p) :
            if not lower_tail :
                return size - i
            else :
                return i
        else :
            sum += p_binom
    if not lower_tail :
        return 0
    else :
        return size
    #
#

# binomial quantile calculation
Dict_qbnom = dict()

# Note : Error will raise if X is too large, for example: >1000
def GetQbnom(p_value, X, prob, low_tail=False) :
    if p_value >=1 :
        return 0
    elif p_value == 0 :
        return X+1
    #
    if X > 1000 :
        X = 1000
    if X not in Dict_qbnom :
        Dict_qbnom[X] = qbinom(p_value, X, prob, low_tail)
    return Dict_qbnom[X]
#

def MatMult(a,b):
    zip_a = zip(*a)
    return [[sum(ele_b*ele_a for ele_b, ele_a in zip(row_b, col_a)) for col_a in zip_a] for row_b in b]
# Example :
#x = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]
#  1  4  7  10
#  2  5  8  11
#  3  6  9  12
#y = [[1,2],[1,2],[3,4]]

#Another way to multiply the vectors
#import numpy as np # I want to check my solution with numpy

#mx = np.matrix([[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
#my = np.matrix([[1,2],[1,2],[3,4]])

"""
# ATCGmap Format Example
chr1    C       4654    CG      CG      0       0       20      0       0       0       0       0       0  1.0
chr1    G       4655    CG      CG      0       0       0       20      0       0       0       0       0  na
chr1    G       4656    CHG     CC      0       0       0       20      0       0       0       0       0  na
chr1    G       4657    CHH     CC      0       0       0       20      0       0       0       0       0  na
chr1    C       4658    CHH     CC      0       19      1       0       0       0       0       0       0  0.05
chr1    C       4659    CHH     CC      0       20      0       0       0       0       0       0       0  0.0
chr1    C       4660    CHH     CC      0       20      0       0       0       0       0       0       0  0.0
chr1    C       4661    CHH     CT      0       19      1       0       0       0       0       0       0  0.05
chr1    T       4662    --      --      0       20      0       0       0       0       0       0       0  na
chr1    C       4663    CHH     CA      0       20      0       0       0       0       0       0       0  0.0
"""



Dict      = {"A":0, "T":1, "C":2, "G":3, "N":4, "T/C":5, "A/G":6}
RevDict   = ["A", "T", "C", "G", "N", "T/C", "A/G"]
VagueDict = {"A":6, "T":5, "C":5, "G":6, "N":4}



AllCases  = {
    "A" : ["A"],
    "T" : ["T"],
    "C" : ["C"],
    "G" : ["G"],
    "T/C" : ["T", "C", "T,C"],
    "A/G" : ["A", "G", "A,G"],
    "A,T" : ["A,T", "T,A"],
    "A,C" : ["A,C", "C,A"],
    "A,G" : ["A,G", "G,A"],
    "T,A" : ["A,T", "T,A"],
    "T,C" : ["T,C", "C,T"],
    "T,G" : ["T,G", "G,T"],
    "C,A" : ["A,C", "C,A"],
    "C,T" : ["C,T", "T,C"],
    "C,G" : ["C,G", "G,C"],
    "G,A" : ["A,G", "G,A"],
    "G,T" : ["G,T", "T,G"],
    "G,C" : ["G,C", "C,G"],
    "A,T/C" : ["A,T", "A,C", "T,A", "C,A"],
    "T,T/C" : ["T", "T,C", "C,T"],
    "C,T/C" : ["C", "T,C", "C,T"],
    "G,T/C" : ["G,T", "G,C", "C,G", "T,G"],
    "A,A/G" : ["A", "A,G", "G,A"],
    "T,A/G" : ["T,A", "T,G", "A,T", "G,T"],
    "C,A/G" : ["C,A", "C,G", "A,C", "G,C"],
    "G,A/G" : ["G", "G,A", "A,G"]
}

UniqCases  = {
    "A" : ["A"],
    "T" : ["T"],
    "C" : ["C"],
    "G" : ["G"],
    "T/C" : ["T", "C", "T,C"],
    "A/G" : ["A", "G", "A,G"],
    "A,T" : ["A,T"],
    "A,C" : ["A,C"],
    "A,G" : ["A,G"],
    "T,A" : ["A,T"],
    "T,C" : ["T,C"],
    "T,G" : ["T,G"],
    "C,A" : ["A,C"],
    "C,T" : ["T,C"],
    "C,G" : ["C,G"],
    "G,A" : ["A,G"],
    "G,T" : ["T,G"],
    "G,C" : ["C,G"],
    "A,T/C" : ["A,T", "A,C"],
    "T,T/C" : ["T", "T,C"],
    "C,T/C" : ["C", "T,C"],
    "G,T/C" : ["C,G", "T,G"],
    "A,A/G" : ["A", "A,G"],
    "T,A/G" : ["T,G", "A,T"],
    "C,A/G" : ["A,C", "C,G"],
    "G,A/G" : ["G", "A,G"]
}

#
# ==========================================================

global binom_options

def PredictNT_binom ( W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G, nuc="N"):
    global binom_options
    p_value = binom_options['pv']
    error_rate = binom_options['er']
    least_cov = binom_options['cov']
    #
    VW = [W_A, W_T, W_C, W_G]
    VC = [C_A, C_T, C_C, C_G]
    # Step 1: The minimum read counts
    SUM_W = sum(VW)
    SUM_C = sum(VC)
    WD_lst = []
    CD_lst = []
    D_lst = []
    # Step 2: Get the LIMIT_W and LIMIT_C
    LIMIT_W = GetQbnom(p_value, SUM_W, error_rate, False)
    if SUM_W >= least_cov :
        WVS = [VW[0], VW[1], VW[2]+VW[1], VW[3]]
        WD_lst = [i for i in range(len(WVS)) if WVS[i] >= LIMIT_W]
    LIMIT_C = GetQbnom(p_value, SUM_C, error_rate, False)
    if SUM_C >= least_cov :
        CVS = [VC[0], VC[1], VC[2], VC[3]+VC[0]]
        CD_lst = [i for i in range(len(CVS)) if CVS[i] >= LIMIT_C]
    # Step 3 : Get the intersection
    if WD_lst != [] :
        if CD_lst != [] :
            D_lst = [ i for i in WD_lst if i in CD_lst ]
        else :
            D_lst = WD_lst
            if 1 in D_lst :
                D_lst.remove(1) # remove T
                D_lst.append(5) # add T/C
            if 2 in D_lst and W_C < LIMIT_W :
                D_lst.remove(2)
            if (2 in D_lst) and (5 in D_lst) and len(D_lst)>2 :
                D_lst.remove(5)
            #
        #
    else :
        if CD_lst != [] :
            D_lst = CD_lst
            if 0 in D_lst :
                D_lst.remove(0) # remove A
                D_lst.append(6) # add A/G
            if 3 in D_lst and C_G < LIMIT_C :
                D_lst.remove(3)
            if (3 in D_lst) and (4 in D_lst) and len(D_lst)>2 :
                D_lst.remove(6)
            #
        else :
            SUM_ALL = SUM_W + SUM_C
            ALLVS = [VW[0], VC[1], VW[2] + VC[2], VW[3] + VC[3], 0, VW[1], VC[0]]
            LIMIT_ALL = GetQbnom(p_value, SUM_ALL, error_rate, False)
            ALL_D_lst = [i for i in range(len(ALLVS)) if ALLVS[i] >= LIMIT_ALL]
            SingleNuc = [i for i in ALL_D_lst if i < 4]
            WildcardNuc = [i for i in ALL_D_lst if i > 4]
            if len(SingleNuc) > 1:
                D_lst = SingleNuc
            elif len(SingleNuc) == 1:
                if WildcardNuc == []:
                    D_lst = SingleNuc
                else:
                    if 5 in WildcardNuc:  # T/C
                        CVS = [VC[0], VC[1], VC[2], VC[0] + VC[3]]
                        CD_lst = [i for i in range(len(CVS)) if CVS[i] >= LIMIT_C]
                        if CD_lst != []:
                            WildcardNuc.remove(5)
                    if 6 in WildcardNuc:  # A/G
                        WVS = [VW[0], VW[1], VW[1] + VW[2], VW[3]]
                        WD_lst = [i for i in range(len(WVS)) if WVS[i] >= LIMIT_W]
                        if WD_lst != []:
                            WildcardNuc.remove(6)
                    D_lst = SingleNuc + WildcardNuc
                    #
                #
            else:  # SingleNuc is empty
                if len(WildcardNuc) == 1:
                    D_lst = WildcardNuc
                else:
                    D_lst = []
                    #
                #
            #
        #
    #
    #
    nuc_prob = 0
    if D_lst == [] :
        return [nuc, 1-p_value, 1-p_value]
    else :
        GENO = ",".join([RevDict[i] for i in D_lst])
        if GENO not in AllCases:
            sys.stderr.write("[Warning] %s not found in genotypes\n" % GENO)
            return [nuc, 1-p_value, 1-p_value]
        #
        if nuc in UniqCases[GENO]:
            nuc_prob = (1-p_value)/len(UniqCases[GENO])
        else :
            nuc_prob = p_value/(10-len(UniqCases[GENO]))
        #
        return [GENO, 1-p_value, nuc_prob]
    #
#


# ===================================
# for the bayesian mode
# import math

ClearSet = ["AA", "AT", "AC", "AG",
                  "TT", "TC", "TG",
                        "CC", "CG",
                              "GG"]
#
ClearSetIndex = {"AA": 0, "AT": 1, "AC": 2, "AG": 3,
                          "TT": 4, "TC": 5, "TG": 6,
                                   "CC": 7, "CG": 8,
                                            "GG": 9}
#
#   IUPAC code:
#       Y = C or T , Y = T/C
#       R = A or G , R = A/G
#
VagueSet = [ "Y",  "R",
            "AY", "TY", "CY", "GY",
            "AR", "TR", "CR", "GR"]
#
VagueSetIndex = { "Y": 0,  "R": 1,
                 "AY": 2, "TY": 3, "CY": 4, "GY": 5,
                 "AR": 6, "TR": 7, "CR": 8, "GR": 9}
#
# ========== A  T  C  G
PriorProb = [1, 2, 2, 2,  # A
                1, 2, 2,  # T
                   1, 2,  # C
                      1]  # G
# ==========    A  T  C  G
LogPriorProb = [0, 1, 1, 1,  # A
                   0, 1, 1,  # T
                      0, 1,  # C
                         0]  # G
#

import math

def CovToPv (Cov) : # for --bayes-dynamicP
    if Cov < 22 :
        return math.exp((Cov-22)/10)
    #
    return 0.99
#

global bayes_options

def PredictNT_bayes ( W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G, nuc="N" ):
    #global p_value
    #global error_rate
    #global dynamicP
    global bayes_options
    p_value = bayes_options['pv']
    error_rate = bayes_options['er']
    dynamicP = bayes_options['dynamicP']
    #
    #Q = 0.01;
    Q = error_rate
    log_Q = math.log(Q) # p(A|TT), p(C|TT), p(G|TT) # error rate
    P = 1-3*Q;      log_P = math.log(P) # p(T|TT)
    q = Q;          log_q = math.log(q) # p(C|AT), p(G|AT) # error rate
    p = (P+Q)/2;    log_p = math.log(p) # p(A|AT)=[p(A|A)+p(A|T)]/2, p(T|AT)=[p(T|A)+p(T|T)]/2
    QQ = Q*2;       log_QQ = math.log(QQ)
    pq = p+q;       log_pq = math.log(pq)
    PQ = P+Q;       log_PQ = math.log(PQ)
    pp = p*2;       log_pp = math.log(pp)
    qq = q*2;       log_qq = math.log(qq)
    #
    # Matrix for Tw
    # ----+----------------------------------
    #  Tw |    A   |   T   |   C    |   G    |
    # ----+----------------------------------|
    #  A  |    2Q  |  p+q  |  p+q   |   2q   |
    #  T  |    -   |  P+Q  |   2p   |   p+q  |
    #  C  |    -   |   -   |  P+Q   |   p+q  |
    #  G  |    -   |   -   |   -    |   2Q   |
    # ----+----------------------------------
    #
    Factor_Tw = [ log_QQ, log_pq, log_pq, log_qq,
                          log_PQ, log_pp, log_pq,
                                  log_PQ, log_pq,
                                          log_QQ]
    #
    # Matrix for Ac
    # ----+----------------------------------
    #  Ac |    A   |   T   |   C    |   G    |
    # ----+----------------------------------|
    #  A  |   P+Q  |  p+q  |  p+q   |   2p   |
    #  T  |    -   |  2Q   |   2q   |   p+q  |
    #  C  |    -   |   -   |   2Q   |   p+q  |
    #  G  |    -   |   -   |   -    |   P+Q  |
    # ----+----------------------------------
    #
    Factor_Ac = [ log_PQ, log_pq, log_pq, log_pp,
                          log_QQ, log_qq, log_pq,
                                  log_QQ, log_pq,
                                          log_PQ]
    #
    ##                 W_A,    W_C,    W_G,    C_T,    C_C,    C_G
    #OtherFactor = [ log_P,  log_Q,  log_Q,  log_Q,  log_Q,  log_Q, # "AA" [0]
    #                log_p,  log_q,  log_q,  log_p,  log_q,  log_q, # "AT" [1]
    #                log_p,  log_p,  log_q,  log_q,  log_p,  log_q, # "AC" [2]
    #                log_p,  log_q,  log_p,  log_q,  log_p,  log_p, # "AG" [3]
    #                #
    #                log_Q,  log_Q,  log_Q,  log_P,  log_Q,  log_Q, # "TT" [4]
    #                log_q,  log_p,  log_q,  log_p,  log_q,  log_q, # "TC" [5]
    #                log_q,  log_q,  log_p,  log_p,  log_q,  log_p, # "TG" [6]
    #                #
    #                log_Q,  log_P,  log_Q,  log_Q,  log_P,  log_Q, # "CC" [7]
    #                log_q,  log_p,  log_p,  log_q,  log_p,  log_p, # "CG" [8]
    #                #
    #                log_Q,  log_Q,  log_P,  log_Q,  log_Q,  log_P, # "GG" [9]
    #              ]
    #               AA,     AT,     AC,     AG,     TT,     TC,     TG,     CC,     CG,     GG
    Factor_Aw = [   log_P,  log_p,  log_p,  log_p,  log_Q,  log_q,  log_q,  log_Q,  log_q,  log_Q ]
    Factor_Cw = [   log_Q,  log_q,  log_p,  log_q,  log_Q,  log_p,  log_q,  log_P,  log_p,  log_Q ]
    Factor_Gw = [   log_Q,  log_q,  log_q,  log_p,  log_Q,  log_q,  log_p,  log_Q,  log_p,  log_P ]
    Factor_Tc = [   log_Q,  log_p,  log_q,  log_q,  log_P,  log_p,  log_p,  log_Q,  log_q,  log_Q ]
    Factor_Cc = [   log_Q,  log_q,  log_p,  log_q,  log_Q,  log_p,  log_q,  log_P,  log_p,  log_Q ]
    Factor_Gc = [   log_Q,  log_q,  log_q,  log_p,  log_Q,  log_q,  log_p,  log_Q,  log_p,  log_P ]
    #
    #
    GenotypeCode = { "AA":"A", "AT":"A,T", "AC":"A,C", "AG":"A,G", "TT":"T", "TC":"T,C", "TG":"T,G",
                     "CC":"C", "CG":"C,G", "GG":"G",
                     "Y":"T/C",  "R":"A/G",
                     "AY":"A,T/C", "CY":"C,T/C", "GY":"G,T/C", "TY":"T,T/C",
                     "AR":"A,A/G", "CR":"C,A/G", "GR":"G,A/G", "TR":"T,A/G", "-":"-"
                   }
    #
    # Bayesian
    #        |   A    T    C    G
    # -----------------------------
    # Watson |  W_A, W_T, W_C, W_G,
    # Crick  |  C_A, C_T, C_C, C_G,
    #
    #[W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G] = [int(i) for i in [ W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G]]
    COV = sum([W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G])
    log_pstP    = [0] * 10
    for i in xrange(10):
        log_pstP[i]   = LogPriorProb[i]  + \
                        W_A*Factor_Aw[i] + W_T*Factor_Tw[i] + W_C*Factor_Cw[i] + W_G*Factor_Gw[i] + \
                        C_A*Factor_Ac[i] + C_T*Factor_Tc[i] + C_C*Factor_Cc[i] + C_G*Factor_Gc[i]
    #
    # ==================================================
    sum_log_pstP = sum( 2**i for i in log_pstP )
    if sum_log_pstP == 0 :
        sys.stderr.write("[warning] line was skipped.\n")
        return [nuc, 1, 0]
    #
    pstP = [2**i/sum_log_pstP for i in log_pstP ]
    VagueSet_pstP = [ pstP[4]+pstP[5]+pstP[7], # Y = TT | TC | CC
                      pstP[0]+pstP[3]+pstP[9], # R = AA | AG | GG
                      pstP[1]+pstP[2], # AY = AT | AC
                      pstP[4]+pstP[5], # TY = TT | TC
                      pstP[5]+pstP[7], # CY = TC | CC
                      pstP[6]+pstP[8], # GY = TG | CG
                      pstP[0]+pstP[3], # AR = AA | AG
                      pstP[1]+pstP[6], # TR = AT | TG
                      pstP[2]+pstP[8], # CR = AC | CG
                      pstP[3]+pstP[9]  # GR = AG | GG
                    ]
    #for i in VagueSet_pstP:
    #    print i
    #
    VagueFactor = [0.93] * 2 + [0.975] * 8
    V_pstP = [i*j for i,j in zip(VagueSet_pstP, VagueFactor)]
    #
    [Prob, geno] = max(zip(pstP, ClearSet) + zip(V_pstP, VagueSet) )
    if nuc in ["A", "T", "C", "G"] :
        prob_nuc = pstP[ClearSetIndex[nuc*2]]
    else :
        prob_nuc = 1.0
    #
    if dynamicP :
        p_value = CovToPv(COV)
    #
    if Prob >= (1-p_value) :
        return [GenotypeCode[geno], Prob, prob_nuc]
    else :
        return [nuc, prob_nuc, prob_nuc]
    #
    # Prob: probability for predicted nuc
    # prob_nuc: probability for NUC
# end

"""
PredictNT_bayes( 3, 3, 0, 0,    0, 0, 0, 0)
"""

import gzip


def VCF_line (CHR, POS, REF, GN, DP, prob_pre, prob_nuc, FILTER = "PASS"):
    CHR = CHR.replace("chr", "").replace("CHR", "").replace("Chr", "")
    ID = "."
    vague = False
    genotype = GN.replace('T/C', 'R').replace('A/G', 'Y')
    ALT_lst = genotype.split(",")
    if 'R' in ALT_lst :
        vague = True
        GU='T/C'
    elif 'Y' in ALT_lst :
        vague = True
        GU='A/G'
    #
    ALT_lst.remove(REF) if REF in ALT_lst else None
    #print ALT_lst
    ALT = ",".join(ALT_lst)
    if ALT == "" :
        ALT = "."
    #
    #print ALTerror_rate
    if prob_pre <1:
        QUAL = int(- math.log10(1 - prob_pre) * 10)
    else :
        QUAL = 99
    #
    if REF in AllCases[GN]:  # no mutation
        GQ = 99  # log10 (0) be very large
    else:
        if (1 - prob_nuc - prob_pre) > 0 :
            GQ = int(- math.log10((1 - prob_nuc - prob_pre) / (1 - prob_nuc)) * 10)
        else :
            GQ = 99
        #
    #
    INFO = ":".join(["NS=1", "DP=%s" % DP])
    if vague :
        INFO = INFO + ":GU=" + GU
    #
    if "," not in genotype :
        GT = "1/1"
    else:
        if "," in ALT :
            GT = "1/2"
        else :
            GT = "0/1"
        #
    #
    # QUAL: quality: Phred-scaled quality score for the assertion made in ALT. i.e.
    #       -10log10 prob(call in ALT is wrong).
    #      If ALT is '.' (no variant) then this is -10log10 prob(variant)
    #      If ALT is not '.' (variant) then the is -10log10 prob(no variant)
    # GQ: conditional genotype quality, encoded as a phred quality -10log10 P(genotype call is wrong,
    #       conditioned on the sites' being variant) (Integer)
    #
    if '/' in GN and REF in AllCases[GN] :
        FILTER = "Vague"
    GENOME = ":".join(["%s"%GT, "%d"%GQ, "%d"%DP])
    return "\t".join([CHR, POS, ID, REF, ALT, "%s"%QUAL, FILTER, INFO, "GT:GQ:DP", GENOME])
#

def SNVFromATCGmap (infile, vcffile, show_all, mode="binom"):
    #
    # global p_value
    # global error_rate
    #
    # Decide the source of input
    if infile :
        if infile.endswith(".gz") :
            IN = gzip.open(infile, 'rb')
        else :
            IN = open(infile, 'r')
        #
    else :
        IN = sys.stdin
    #
    # Decide the mode for calling SNP
    if mode == "binom" :
        SNPfunc = PredictNT_binom
        global binom_options
        sys.stderr.write("# BinomWC mode\n")
        sys.stderr.write("# options: p value = %f\n" % binom_options["pv"])
        sys.stderr.write("# options: error rate = %f\n" % binom_options["er"])
        sys.stderr.write("# options: checkpoint coveraeg = %d\n" % binom_options["cov"])
    elif mode =="bayes" :
        SNPfunc = PredictNT_bayes
        global bayes_options
        sys.stderr.write("# BayesWC mode\n")
        if bayes_options["dynamicP"] :
            sys.stderr.write("# options: use dynamic p-values for different coverages\n" )
        else :
            sys.stderr.write("# options: p value = %f\n" % bayes_options["pv"])
        #
        sys.stderr.write("# options: error rate = %f\n" % bayes_options["er"])
    else :
        sys.stderr.write("[Error] Wrong mode specified.\n")
        return None
    #
    #
    # VCF format
    VCF = None
    if vcffile :
        if vcffile.endswith(".gz") :
            VCF = gzip.open(vcffile, 'wb')
        else :
            VCF = open(vcffile, 'w')
        #
        VCF.write("##fileformat=VCFv4.2\n")
        from datetime import datetime
        DateStr=datetime.now().strftime('%Y%m%d')
        VCF.write("##fileDate=%s\n" % DateStr)
        VCF.write("##source=CGmapTools\n")
        VCF.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
        VCF.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        VCF.write('##INFO=<ID=GU,Number=.,Type=String,Description="Genotype is Unidentified:R=A/G, Y=T/C">\n')
        VCF.write('##FILTER=<ID=q10,Description="Quality below 10">\n')
        VCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        VCF.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
        VCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        VCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n")
    #
    print "\t".join(["#chr", "nuc", "pos", "ATCG_watson", "ATCG_crick", "predicted_nuc", "p_value"] )
    #
    for line in IN :
        if not line.strip() :
            continue
        else :
            try :
                chr, nuc, pos, pattern, dinuc, W_A, W_T, W_C, W_G, W_N, C_A, C_T, C_C, C_G, C_N, methyl = line.strip().split()
            except ValueError :
                sys.stderr.write("\n[Error]:\n\t File [ %s ] may have wrong number of columns.\n" % infile)
                exit(-1)
            # end_of_try
            [W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G] = [int(i) for i in [W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G]]
            [NUC, prob_pre, prob_nuc] = SNPfunc( W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G, nuc)
            COV = sum([W_A, W_T, W_C, W_G, C_A, C_T, C_C, C_G])
            #
            if show_all :
                print "\t".join([chr, nuc, pos,
                                 "%d, %d, %d, %d"%(W_A, W_T, W_C, W_G), "%d, %d, %d, %d"%(C_A, C_T, C_C, C_G),
                                 NUC, "%.2e" % (1-prob_pre) ] )
                if vcffile:
                    VCF.write(VCF_line(chr, pos, nuc, NUC, COV, prob_nuc, prob_nuc) + "\n" )
            elif NUC != "-" :
                #
                if (NUC in AllCases) :
                    if (nuc not in AllCases[NUC]) :
                        print "\t".join([chr, nuc, pos,
                                        "%d,%d,%d,%d"%(W_A, W_T, W_C, W_G), "%d,%d,%d,%d"%(C_A, C_T, C_C, C_G),
                                         NUC, "%.1e" % (1-prob_pre) ] )
                        if vcffile :
                            VCF.write(VCF_line(chr, pos, nuc, NUC, COV, prob_nuc, prob_nuc) + "\n" )
                    #
                else :
                    sys.stderr.write( "[Warning] %s not found in genotypes\n" %  NUC )
                #
            #
        #
    # for
    if infile :
        IN.close()
    #
#

from optparse import OptionParser
#from textwrap import dedent

# ===========================================
def main():
    usage = "Usage: cgmaptools snv [-i <ATCGmap>] [-o <output> -v <VCF>]\n" \
            "      (aka SNVFromATCGmap)\n" \
            "Description: Predict the SNV from ATCGmap file.\n" \
            "Contact:     Guo, Weilong; guoweilong@126.com\n" \
            "Last update: 2016-12-07\n" \
            "Output format example:\n" \
            "   #chr  nuc  pos    ATCG_watson  ATCG_crick  predicted_nuc  p_value\n" \
            "   chr1  G    4752   17,0,0,69    0,0,0,0     A,G            9.3e-07\n" \
            "   chr1  A    4770   40,0,0,29    0,0,0,0     A,G            0.0e+00\n" \
            "   chr1  T    8454   0,39,0,0     0,0,0,0     T/C            1.00e-01\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile", help="ATCGmap format, STDIN if not specified",
                      metavar="FILE", default=None)
    parser.add_option("-v", "--vcf", dest="vcffile", help="VCF format file for output",
                      metavar="FILE", default=None)
    parser.add_option("-a", "--all_nt", action="store_true", dest="showall", default = False,
                      help = 'Show all sites with enough coverage (-l). Only show SNP sites if not specified.')
    parser.add_option("-o", dest="outfile", default=None, help="STDOUT if not specified")
    parser.add_option("-m", "--mode", dest="mode",
                      help="Mode for calling SNP [Default: %default]\
                            binom: binomial,  separate strands\
                            bayes: bayesian mode",
                      default="binom")
    # options for BayesWC mode
    parser.add_option("--bayes-e", dest="bayes_er",
                      help="(BayesWC mode) Error rate for calling a nucleotide [Default: %default]",
                      type="float", default=0.05)
    parser.add_option("--bayes-p", dest="bayes_pv",
                      help="(BayesWC mode) P value as cut-off [Default: %default]",
                      type="float", default=0.001)
    parser.add_option("--bayes-dynamicP", action="store_true", dest="bayes_dynamicP", default = False,
                      help = "(BayesWC mode) Use dynamic p-value for different coverages install of "
                             "specific p-value. (Recomended)\
                              \"--bayes-p\" will be ignored if \"--bayes-dynamicP\" is specified.")
    # options for BinomWC mode
    parser.add_option("--binom-e", dest="binom_er",
                      help="(BinomWC mode) Error rate for calling a nucleotide [Default: %default]",
                      type="float", default=0.05)
    parser.add_option("--binom-p", dest="binom_pv",
                      help="(BinomWC mode) P value as cut-off [Default: %default]",
                      type="float", default=0.01)
    parser.add_option("--binom-cov", dest="binom_cov",
                      help="(BinomWC mode) The coverage checkpoint [Default: %default]",
                      type="int", default=10)
    (options, args) = parser.parse_args()
    #

    if (options.outfile is not None) :
        sys.stdout = open(options.outfile, 'w')
    #
    #global p_value
    #global error_rate
    #global least_cov
    #global dynamicP
    #
    #
    # ======================
    global bayes_options
    bayes_options = {"dynamicP": True,
                     "pv": 0.001,
                     "er": 0.05}
    bayes_options["pv"] = float(options.bayes_pv)
    bayes_options["er"] = float(options.bayes_er)
    bayes_options["dynamicP"] = options.bayes_dynamicP
    #
    # ======================
    global binom_options
    binom_options = {"pv": 0.001,
                     "er": 0.05,
                     'cov': 5}
    binom_options['pv'] = float(options.binom_pv)
    binom_options['er'] = float(options.binom_er)
    binom_options['cov'] = int(options.binom_cov)
    #
    #p_value = float(options.p_value)
    #error_rate = float(options.error_rate)
    #least_cov = int(options.least_cov)
    #dynamicP = options.dynamicP
    #
    SNVFromATCGmap( options.infile, options.vcffile, options.showall,
                    options.mode)
    #

# ===========================================
if __name__ == "__main__":
    main()

