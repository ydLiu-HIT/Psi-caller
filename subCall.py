#!/usr/bin/env python
# coding=utf-8

###
#  get genotype likelihood by realignment mapping score
###

import os
import sys
import gc
import logging
import math
import re
#import numpy as np
from time import time
from utils import base_to_ACGT as base2ACGT 
from check import *
from ksw_module import ksw_aligner
from ksw2_module import ksw2_aligner
from math_func import *
from Canclass import *
from getHCVariant import callLCVariant

majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
_MIN_MQ = 1e-10
_MAX_LIKELIHOODSCORE = 1e-250
_MAX_PROB = 1.0

EXTFLANK = 5
AMPLIFIER = 10.0

CDEBUG = False
DEBUG = False 

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def evc_base_from(base):
    return base if base == "N" else base2ACGT[base]

#def haveMatch(pos, typ, IndelCan, shift, vbase):
#    mindis = shift; minid = 0; sig = False
#    for i in range(len(IndelCan)):
#        r = int(IndelCan[i][0])
#        t = IndelCan[i][1]
#        if abs(pos - r) < mindis:
#            sig = True
#            break
#        #if mindis > 0:
#        #    if abs(pos - r) < mindis and typ == t:
#        #    #if abs(pos - r) < mindis:
#        #        mindis = abs(pos-r); minid = i; sig = True
#        #else:
#        #    if pos == r and typ == t and vbase == IndelCan[i][4]:
#        #        minid = i; sig = True
#
#    return sig, minid

def haveMatch(pos, typ, IndelCan, shift, vbase):
    mindis = shift; minid = 0; sig = False
    for i in range(len(IndelCan)):
        r = int(IndelCan[i][0])
        t = IndelCan[i][1]
        if mindis > 0:
            if abs(pos - r) < mindis and typ == t:
                mindis = abs(pos-r); minid = i; sig = True
        else:
            if pos == r and typ == t and vbase == IndelCan[i][4]:
                minid = i; sig = True

    return sig, minid


def makeCandidate(chrName, alignment, RNAME, REF, rS, SEQ, shift, source, multiCan, max_merge_dis):
    quePos = 0; TS = int(RNAME[2]); refPos = TS
    IndelCan = []; Candidates = []
    for i in RNAME[4:]:
        IndelCan.append(i.split('|'))
    leftmost = int(IndelCan[0][0]) - max_merge_dis
    rightmost = int(IndelCan[-1][0]) + max_merge_dis
    LL = len(IndelCan); idx = 0
    CIGAR = alignment[1]
    op_l = 0
    for m in str(CIGAR)[:-1]:
        if m.isdigit():
            op_l = op_l * 10 + int(m)
            continue
        if m == "S":
            quePos += op_l
        elif m == "M" or m == "=": # M, X, =
            refPos += op_l
            quePos += op_l
        elif m == "X":
            for _ in range(op_l):
                refPos += 1; quePos += 1;
                rseq = evc_base_from(REF[refPos-rS]); qseq = evc_base_from(SEQ[quePos-1])
                if rseq == qseq: continue
                if multiCan:
                    if refPos > leftmost and refPos < rightmost:
                        cand = [chrName,refPos,"S",rseq,qseq,source,alignment[0],0,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                else:
                    sig, idx = haveMatch(refPos, "S", IndelCan, 2, qseq)
                    if sig == True:
                        cand = [chrName,refPos,"S",rseq,qseq,source,alignment[0],0,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                        del IndelCan[idx]
        elif m == "I":
            if op_l < 55:
                rseq = evc_base_from(REF[refPos-rS]); qseq = SEQ[quePos-1:quePos+op_l]
                if multiCan:
                    if refPos > leftmost and refPos < rightmost:
                        cand = [chrName,refPos,"INS",rseq,qseq,source,alignment[0],op_l,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                else:
                    sig, idx = haveMatch(refPos, "I", IndelCan, shift, SEQ[quePos:quePos+op_l])
                    if sig == True:
                        cand = [chrName,refPos,"INS",rseq,qseq,source,alignment[0],op_l,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                        del IndelCan[idx]
            quePos += op_l
        elif m == "D":
            if op_l < 55:
                rseq = REF[refPos-rS:refPos+op_l+1-rS]; qseq = evc_base_from(SEQ[quePos-1])
                if multiCan:
                    if refPos > leftmost and refPos < rightmost:
                        cand = [chrName,refPos,"DEL",rseq,qseq,source,alignment[0],op_l,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                else:
                    sig, idx = haveMatch(refPos, "D", IndelCan, shift, REF[refPos-rS:refPos+op_l-rS])
                    if sig == True:
                        cand = [chrName,refPos,"DEL",rseq,qseq,source,alignment[0],op_l,int(IndelCan[idx][2]),int(IndelCan[idx][3]), quePos, refPos-TS, SEQ]
                        Candidates.append(cand)
                        del IndelCan[idx]
            refPos += op_l
        op_l = 0
    
    return Candidates 

def re_calc_MQ_locally(cigars, qual, ts, vl, extflank):
    start = ts - extflank
    end = ts + extflank + vl
    rpos = 0; qpos = 0; prob = 0
    if CDEBUG: print("cigar:", cigars)
    advance = 0
    for op in cigars:
        if op.isdigit():
            advance = advance * 10 + int(op)
            continue
        if op in ("=", "X", "D"):
            tp = advance + rpos - start
            rpos += advance; qpos += 0 if op == "D" else advance
            if tp > 0:
                tp = tp if rpos < end else tp+end-rpos
                qpos = qpos if rpos < end else qpos if op == "D" else qpos+end-rpos
                if op == "=":
                    if CDEBUG: print("=", qpos, tp, end=' ')
                    for i in range(qpos-tp, qpos):
                        Q = ord(qual[i])-33 
                        Q = Q if Q > 0 else 1
                        if CDEBUG: print(Q, end=' ')
                        prob += math.log10(1 - math.pow(10,-Q/10.0))
                    if CDEBUG: print('\n')
                elif op == "X":
                    if CDEBUG: print("X", qpos, tp, end=' ')
                    for i in range(qpos-tp, qpos):
                        Q = ord(qual[i])-33 
                        Q = Q if Q > 0 else 1
                        if CDEBUG: print(Q,)
                        prob += (-Q/10.0)
                    if CDEBUG: print('\n')
                elif op == "D":
                    if CDEBUG: print("D", qpos, tp, end=' ')
                    Q_max = 0
                    s = qpos-2 if qpos-2>0 else 0
                    e = qpos+3 if qpos+3<len(qual) else len(qual)
                    for i in range(s, e):
                        Q = ord(qual[i])-33
                        Q = Q if Q > 0 else 1
                        Q_max = Q if Q > Q_max else Q_max
                    prob += (-Q_max / 10.0) * tp
                start = rpos
            if rpos >= end:
                break
        elif op == "I":
            if rpos < start:
                qpos += advance
                advance = 0
                continue
            if CDEBUG: print("I", qpos, advance,)
            for i in range(qpos, qpos + advance):
                Q = ord(qual[i])-33
                Q = Q if Q > 0 else 1
                if CDEBUG: print(Q, end=' ')
                prob += (-Q/10.0)
            qpos += advance
            if CDEBUG: print("\n")
        advance = 0
    if CDEBUG: print("prob=", prob) 
    return prob

def getExactMatch(Cigars, target, query, start):
    Ts = start; Qs = 0
    newCigar = ''
    advance = 0
    for op in Cigars:
        if op.isdigit():
            advance = advance * 10 + int(op)
            continue
        if op == "D":
            Ts += advance
            newCigar += (str(advance)+op)
        elif op == "I" or op =="S":
            Qs += advance
            newCigar += (str(advance)+op)
        elif op == "M":
            tM = 0; tX = 0
            for i in range(advance):
                if target[Ts+i] == query[Qs+i]:
                    if tX > 0: newCigar += (str(tX)+"X"); tX = 0
                    tM += 1
                else:
                    if tM > 0: newCigar += (str(tM)+"="); tM = 0
                    tX += 1

            newCigar += (str(tM)+"=") if tM > 0 else (str(tX)+"X")
            
            Ts += advance; Qs += advance
        advance = 0
    return newCigar


def realignWithQuality(QUERY, QUAL, Haplot, var, Param, amplifier, extflank, x_score):
    if DEBUG:
        print("amplifier = ", amplifier)
    likeliScore = []
    for m in range(len(QUERY)):
        qseq = bytes.decode(QUERY[m])
        qual = QUAL[m]

        rsc = []; od = 0
        for hap in Haplot:
            TS, TL, LH = Param[od]; od += 1;
            alignment = ksw_aligner(qseq, hap, x_score)
            newCigar = getExactMatch(alignment[1], hap, qseq, 0); alignment[1] = newCigar
            prob = re_calc_MQ_locally(alignment[1], qual, TS, TL, extflank); alignment[0] = prob
            alignment[0] = min(max(math.pow(10, prob/amplifier), _MIN_MQ), _MAX_PROB) # method1
            rsc.append(alignment)
        else:
            likeliScore.append(rsc)

    return likeliScore

def realignToHaplotype(QUERY, QUAL, Haplot, var, x_score, Sig): 
    Param = []
    VLen = str(var.varlen).split(',')
    QPos = str(var.quePos).split(',')
    Vtype = var.type.split(',')
    #********************************************
    extflank = 1 if var.type == "S" else EXTFLANK
    #********************************************
    TS = var.refPos; TL = max([int(i) for i in VLen]) if "DEL" in Vtype else 0; LH = extflank*2 + TL
    Param.append((TS, TL, LH))
    for l in range(len(VLen)):
        TS_alt = int(QPos[l]); TL_alt = int(VLen[l]) if Vtype[l] == "INS" else 0; LH_alt = extflank*2 + TL_alt
        Param.append((TS_alt, TL_alt, LH_alt))
    #get local reads
    if DEBUG:
        print(Param)

    #amplifier = 3.5 if var.type == 'S' else AMPLIFIER
    amplifier = 2.0 if var.type == 'S' else AMPLIFIER
    likeliScore = realignWithQuality(QUERY, QUAL, Haplot, var, Param, amplifier, extflank, x_score)

    return likeliScore

def most_likely_genotype(CandidateList, QUERY, QUAL, x_score, Sig):
    for var in CandidateList:
        if Sig:
            print(var.pos, var.type, var.ALT, var.varlen, var.quePos, var.refPos, var.ref, var.alt)
        #get halplotypes
        Haplot = []
        Haplot.append(var.ref)  # ref
        Haplot.extend(var.alt)

        if Sig:
            print("ref:", Haplot[0])
            print("alt:", Haplot[1:])
        
        #Based on alignment score #
        likeliScore = realignToHaplotype(QUERY, QUAL, Haplot, var, x_score, Sig)
        #if DEBUG:
        #    print(np.array(likeliScore))
        GenotypeLikeliHood = getGenotypeLikeliHood(likeliScore, len(Haplot))
        log10_probs = normalize_log10_probs([math.log10(p) for p in GenotypeLikeliHood])
        real_probs = toRealSpace(log10_probs)
        #predictions = round_gls(real_probs)
        var.GT, var.GQ, var.qual = computer_quals(real_probs, len(Haplot))
        
        if Sig:
            print("genotype likelihood:", GenotypeLikeliHood)
            print("log10_prob", log10_probs)
            print("real_prob:", real_probs)
        
def confidence_read(prob):
    pV = [p[0] for p in prob]
    Midx = pV.index(max(pV))
    midx = pV.index(min(pV))
    #if pV[Midx] - pV[midx] < 0.5:
    #    return None
    if pV[Midx] < pV[midx] * 3 and pV[Midx] > 0.99 and pV[midx] < 0.7:  #0.999 get max 
        prob[midx][0] = 0.1

    return prob

def getGenotypeLikeliHood(likeliScore, H):
    HL = sum(x+1 for x in range(H))
    GenotypeLikeHood = [1.0 for i in range(HL)]
    breaksig = False

    for prob in likeliScore:
        prob = confidence_read(prob)

        gt = 0
        for i in range(H):
            for j in range(0, i + 1):
                p = (prob[i][0] + prob[j][0])/2
                if i != j: GenotypeLikeHood[gt] *= p if p < 0.6 else 0.6
                else: GenotypeLikeHood[gt] *= p
                gt += 1

        for t in GenotypeLikeHood:
            if t < _MAX_LIKELIHOODSCORE:
                breaksig = True; break
        if breaksig:
            break
    
    return GenotypeLikeHood

def createVariantCandidate(value, con1, con2, Tc, ref, isFastq, assty, *alt):
    var = variantsCan(value, con1, con2, Tc, assty)
    
    if isFastq:
        for a in alt:
            var.set_alt(a)

        qpos = str(var.quePos).split(',')[0]
        var.set_ref(get_ref_haplotype(alt[0], int(qpos), var.REF, var.ALT.split(',')[0]))
        var.refPos = int(qpos)
    
    return var

def get_ref_haplotype(alt_hap, pos, REF, ALT):
    ref_hap = alt_hap[0:pos-1]
    ref_hap += REF
    ref_hap += alt_hap[pos+len(ALT)-1:]
    
    return ref_hap


def generate_MNP(ins1, ins2, clu_F, clu_R, Tc, refseq, isFastq, assty):
    if ins1[2] == ins2[2]: #same variant type
        if ins2[4] != ins1[4]:
            ins1[2] = ins1[2]+','+ins2[2]
            ins1[4] = ins1[4] + ',' + ins2[4]
            ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
            if ins1[8] + ins2[8] > ins1[9]: ins2[8] = ins1[9] - ins1[8]
            ins1[8] = str(ins1[8]) + ',' + str(ins2[8])
            ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
            var = createVariantCandidate(ins1[:-1], clu_F, clu_R, Tc, refseq, isFastq, assty, ins1[-1], ins2[-1])
        elif ins1[2] == "DEL" and ins1[3] != ins2[3]:
            ins1[2] = ins1[2]+','+ins2[2]
            ins1[3] = ins1[3] if len(ins1[3]) > len(ins2[3]) else ins2[3]
            ins1[4] = ins1[3][0:len(ins1[3])-ins1[7]] + ',' + ins1[3][0:len(ins1[3])-ins2[7]]
            ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
            if ins1[8] + ins2[8] > ins1[9]: ins2[8] = ins1[9] - ins1[8]
            ins1[8] = str(ins1[8]) + ',' + str(ins2[8])
            ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
            var = createVariantCandidate(ins1[:-1], clu_F, clu_R, Tc, refseq, isFastq, assty, ins1[-1], ins2[-1])
        else: # the same variant record  
            var = createVariantCandidate(ins1[:-1], (clu_F+clu_R) if (clu_F+clu_R) < Tc else Tc , 0, Tc, refseq, isFastq, assty, ins1[-1])
            var.supportA = clu_F+clu_R
            var.supportT = Tc
    else: #the same pos, ins and del
        if len(ins1[3]) < len(ins2[3]):
            ins1[3] = ins2[3]
            dl = ins2[3][1:]
            dA = ins1[4]+dl+','+ins2[4]
        else:
            dl = ins1[3][1:]
            dA = ins1[4]+','+ins2[4]+dl

        #ins1[3] = ins1[3] if len(ins1[3]) > len(ins2[3]) else ins2[3]
        #dl = ins1[3][1:]
        #dA = ins1[4]+dl #if len(ins1[4]) > 1 else ins1[4]
        #dA += ','
        #dA += ins2[4]+dl #if len(ins2[4]) > 1 else ins2[4]
        ins1[2] = ins1[2] + ',' + ins2[2]
        ins1[4] = dA
        ins1[7] = str(ins1[7]) + ',' + str(ins2[7])
        ins1[8] = str(ins1[8]) + ',' + str(ins2[8])
        ins1[10] = str(ins1[10]) + ',' + str(ins2[10])
        var = createVariantCandidate(ins1[:-1], clu_F, clu_R, Tc, refseq, isFastq, ins1[-1], ins2[-1])
    
    return var

def mergeCoVar(canList, cluster_ids_n, Tc, refseq, isFastq, assty):
    variantList = list()
    het_n = len(canList)
    if het_n == 2:
        listF, listR = canList
        clu_F, clu_R = cluster_ids_n
        lenF = len(listF); lenR = len(listR)
        idF = 0; idR = 0
        while idF < lenF and idR < lenR:
            varF = listF[idF]; varR = listR[idR]
            if varF[1] < varR[1]:
                var = createVariantCandidate(varF[:-1], clu_F, clu_R, Tc, refseq, isFastq, assty, varF[-1])
                idF += 1
            elif varF[1] > varR[1]:
                var = createVariantCandidate(varR[:-1], clu_R, clu_F, Tc, refseq, isFastq, assty, varR[-1])
                idR += 1
            else:
                var = generate_MNP(varF, varR, clu_F, clu_R, Tc, refseq, isFastq, assty)
                idF += 1; idR += 1

            variantList.append(var)

        while idF < lenF:
            var = createVariantCandidate(listF[idF][:-1], clu_F, clu_R, Tc, refseq, isFastq, assty, listF[idF][-1])
            variantList.append(var)
            idF += 1
        while idR < lenR:
            var = createVariantCandidate(listR[idR][:-1], clu_R, clu_F, Tc, refseq, isFastq, assty, listR[idR][-1])
            variantList.append(var)
            idR += 1
    elif het_n == 1:
        for v in canList[0]:
            var = createVariantCandidate(v[:-1], Tc, 0, Tc, refseq, isFastq, assty, v[-1])
            variantList.append(var)

    return variantList


def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input")
    parser.add_argument("--fout_vcf", type=str, default="detect.vcf",
            help="vcf path for detected variants, default: %(default)s")
    parser.add_argument("--workDir", type=str, default="tempDir",
            help="Temporary working path for Multiple sequence alignment, default:%(default)s")
    parser.add_argument("--canfix", type=str, default="candidate",
            help="Candidate file prefix used in ExtractVariant.py, default: %(default)s")
    parser.add_argument("--shift", type=int, default=5,
            help="The distance bewteen the detect variant and the activate region, default:%(default)d")
    parser.add_argument("--perror_for_snp", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")
    parser.add_argument("--perror_for_indel", type=float, default="0.1",
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    parser.add_argument("--debug", type=bool, default=False,
            help="debugs, default:%(default)s")
    args = parser.parse_args()    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    setup_logging()
    cwd = os.getcwd()
    args.workDir = cwd + '/' + args.workDir + '/'

    global DEBUG
    DEBUG = True if args.debug else False

    main_ctrl(args)

if __name__ == "__main__":
    t_s = time()
    run()
    t_e = time()
    logging.info("Finish the script in %f seconds", t_e - t_s)





