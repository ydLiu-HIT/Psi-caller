#!/usr/bin/env python
# coding=utf-8

import os
import sys
import gc
import logging
import shlex
import argparse
from time import time
from utils import *
from check import *
from getHCVariant import callHCVariant, callLCVariant
from subCall import mergeCoVar, most_likely_genotype, getExactMatch, makeCandidate
from poa_module import POA
from subPoa_module import subPOA
from ksw2_module import ksw2_aligner
from ksw_module import ksw_aligner
from mantaAss_module import mantaAss

is_pypy = '__pypy__' in sys.builtin_module_names

MAXR = 250
is_MSA = 0 

majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})

class candidate_instance():
    def __init__(self, pos, vTyp, vLen, refC, rC, var, loci):
        self.pos = pos
        self.vTyp = vTyp
        self.vLen = vLen
        self.refC = refC
        self.rC = rC
        self.var = var
        self.loci = loci

def PypyGCCollect(signum, frame):
    gc.collect()
    signal.alarm(60)

def can_combine(svT1, svT2, S1, S2):
    if svT1 == 'S' and svT2 == 'S':
        return True
    if svT1 != 'S' and svT2 != 'S':
        return True
    if S2-S1 < 20:
        return True
    return False


def candidate_from_file(fin_can, max_merge_dis, is_range_given, refStart, refEnd):
    is_read_file_from_standard_input = fin_can == "PIPE"
    if is_read_file_from_standard_input:
        candidate_file_path_output = sys.stdin
    else:
        candidate_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (fin_can)))
        candidate_file_path_output = candidate_file_path_process.stdout
    
    confuse_thre = max_merge_dis #20
    windowL = 500
    tmpList = list(); variant = list()
    is_first = True
    Start = 0; pres = 0; varType = None
    for row in candidate_file_path_output:
        line = row.strip().split()
        POS = int(line[1])
        ctg = line[0]
        #if ctg not in majorContigs: continue
        if is_range_given and not (POS > refStart and POS < refEnd): continue

        if is_first:
            is_first = False; pres = POS
            tmpList.append(line)
            continue

        if POS - pres < windowL:
            tmpList.append(line)
            pres = POS
        else:
            pres = POS; Len = len(tmpList); t = 0; S = t
            while t < Len:
                if (t < Len-1) and (int(tmpList[t+1][1])-int(tmpList[t][1]) < confuse_thre):
                    t += 1; continue
                
                if t > S or tmpList[t][-1] != "HC":
                    for c in range(S,t+1):
                        item = tmpList[c]
                        tpos = int(item[1])
                        rC = int(item[2]); refB = item[3]; refC = int(item[4])
                        l = 5
                        while l < len(item)-1:
                            sv = item[l]; cnt = int(item[l+1]); l += 2
                            svT = sv[0]
                            if svT == "I" or svT == "D":
                                svB = sv[1:]
                                candidate = candidate_instance(tpos, svT, len(svB), cnt, rC, svB, item[-1])
                            else:
                                candidate = candidate_instance(tpos, 'S', 0, cnt, rC, svT, item[-1])

                            variant.append(candidate)

                        if variant[-1].pos - variant[0].pos > windowL:
                            yield variant
                            variant = []
                    if len(variant) > 0:
                        yield variant
                        variant = []
                else:
                    yield " ".join(tmpList[t])
                t += 1; S = t
            tmpList = []; tmpList.append(line)

    if len(tmpList) > 0:
        Len = len(tmpList); t = 0; S = t
        while t < Len:
            if (t < Len-1) and (int(tmpList[t+1][1])-int(tmpList[t][1]) < confuse_thre):
                t += 1; continue
            
            if t > S or tmpList[t][-1] != "HC":
                for c in range(S,t+1):
                    item = tmpList[c]
                    tpos = int(item[1])
                    rC = int(item[2]); refB = item[3]; refC = int(item[4])
                    l = 5
                    while l < len(item)-1:
                        sv = item[l]; cnt = int(item[l+1]); l += 2
                        svT = sv[0]
                        if svT == "I" or svT == "D":
                            svB = sv[1:]
                            candidate = candidate_instance(tpos, svT, len(svB), cnt, rC, svB, item[-1])
                        else:
                            candidate = candidate_instance(tpos, 'S', 0, cnt, rC, svT, item[-1])

                        variant.append(candidate)

                    if variant[-1].pos - variant[0].pos > windowL:
                        yield  variant
                        variant = []
                if len(variant) > 0:
                    yield variant
                    variant = []
            else:
                yield " ".join(tmpList[t])
            t += 1; S = t
            
    if not is_read_file_from_standard_input:
        candidate_file_path_output.close()
        candidate_file_path_process.wait()


def get_reference_sequence(ref_path, chrName, chrStart, chrEnd):
    refernce_sequences = []

    if chrStart != None and chrEnd != None:
        region_value = "{}:{}-{}".format(chrName, chrStart, chrEnd)
    else:
        region_value = str(chrName)
    
    samtools_faidx_process = subprocess_popen(
        shlex.split("samtools faidx {} {}".format(ref_path, region_value))
    )

    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence

def getFlanking(candi):
    indL = 0
    ifdel = False
    ifLRP = False
    for ci in candi:
        indL = ci.vLen if indL < ci.vLen else indL
        if ci.vTyp == "D": ifdel = True
        if ci.loci == "LRP": ifLRP = True
    
    if ifdel: 
        return indL, ifLRP
    elif indL == 0: 
        return 0, ifLRP
    else:
        return 1, ifLRP

def checksoft_and_realign(ref, rS, queries, Start, End, indL, Sig):
    for l in queries:
        POS, POSEND, SEQ, CIGAR = l[0:4]
        advance = 0
        leftClip = 0; rightClip = 0; isleft = True

        ## for read with large soft clipping
        for m in str(CIGAR):
            if m.isdigit():
                advance = advance * 10 + int(m)
                continue
            if m == "S":
                if isleft: leftClip = advance 
                else: rightClip = advance
            else:
                isleft = False
            advance = 0
        
        thre = 10; re = False
        sv_sim = 50
        if leftClip > thre and rightClip > thre:
            ps = POS-sv_sim - leftClip; pe = POSEND+sv_sim+rightClip; re = True
        elif leftClip > thre:
            ps = POS-sv_sim-leftClip; pe = POSEND; re = True
        elif rightClip > thre:
            ps = POS; pe = POSEND+sv_sim+rightClip; re = True
    
        if re:
            tseq = str(ref[ps+1-rS:pe+1-rS])
            reali = ksw_aligner(SEQ, tseq, 4)
            if Sig:
                print("preCigar:", CIGAR, POS, POSEND)
                print("query:", SEQ)
                print("target:", tseq)
                print("realign:", reali[1], ps, pe)
            l[0] = ps; l[1] = pe; l[3] = reali[1]


def getLocalSeq(ref, rS, queries, candi, Start, End, useBaseQuality, useAllReads, flank, Sig):
    queryCover = list(); queryPart = list()
    qualCover = list(); qualPart = list()
    target = ''
    cnt = 0
    for l in queries:
        POS, POSEND = l[0:2]
        if useBaseQuality:
            SEQ, CIGAR, QUAL = l[2:]
        else:
            SEQ, CIGAR = l[2:]

        if POS > End or POSEND < Start: continue 

        refStart = POS; readStart = 0
        ExtractStart = readStart; ExtractEnd = 0 
        flagS = False; flagE = True; inner = True; spanS = False; spanE = False; spanning = False; leftClip = True
        if refStart < Start: flagS = True 

        op_l = 0; lastop = 0; lastopl = 0
        for m in str(CIGAR):
            if m.isdigit():
                op_l = op_l * 10 + int(m)
                continue

            if m == "S":
                readStart += op_l
                if leftClip:
                    ExtractStart += op_l
                leftClip = False
            elif m == "M" or m == "=" or m == "X": # M, X, =
                refStart += op_l
                readStart += op_l
                leftClip = False
            elif m == "I": # I
                readStart += op_l
                leftClip = False
            elif m == "D": # D
                refStart += op_l
                leftClip = False
            lastop = m; lastopl = op_l
            op_l = 0
            
            # cal start
            if refStart < Start: continue 
            elif flagS:
                if m == "M" or m == "=" or m == "X":
                    ExtractStart = readStart - (refStart - Start)
                else:
                    ExtractStart = readStart
                flagS = False
                spanS = True
            # cal end
            if refStart >= End:
                if m == "M" or m == "=" or m == "X":
                    ExtractEnd = readStart - (refStart - End)
                else:
                    ExtractEnd = readStart
                flagE = False
                spanE = True
                break
        if flagE:
            if lastop == "S" or lastop == "I":
                ExtractEnd = readStart - lastopl
            else:
                ExtractEnd = readStart

        Lread = ExtractEnd - ExtractStart
        Lref = End - Start

        if spanS and spanE: spanning = True
        
        #if Sig:
        #    print(Start, End, POS, POSEND, ExtractStart, ExtractEnd, Lread, Lref, spanning, spanS, spanE, CIGAR)
        
        if cnt < MAXR and Lread > flank:
            if spanning:
                queryCover.append(bytes(SEQ[ExtractStart: ExtractEnd], 'utf-8'))
                if useBaseQuality:
                    qualCover.append(QUAL[ExtractStart: ExtractEnd])
            else:
                if (spanS and End - refStart < flank / 2) or (spanE and POS - Start < flank / 2) or useAllReads:
                    queryPart.append(bytes(SEQ[ExtractStart: ExtractEnd], 'utf-8'))
                    if useBaseQuality:
                        qualPart.append(QUAL[ExtractStart: ExtractEnd])
                
            cnt += 1
   
    name = "RefSeq_" + "id_" +str(Start) + "_" + str(End) + "_" + "_".join(["%d|%s|%d|%d|%s|%s" %(x.pos,x.vTyp,x.refC,x.rC, x.var, x.loci) for x in candi])
    target = str(ref[Start+1-rS:End+1-rS])

    return queryCover+queryPart, qualCover+qualPart, target, name

def ExtractSeqAroundVariant(fin_bam, chrStart, chrEnd, chrName, minMQ, useBaseQuality):
    if chrStart != None and chrEnd != None:
        region_value = "{}:{}-{}".format(chrName, chrStart, chrEnd)
    else:
        region_value = str(chrName)
    p2 = subprocess_popen(shlex.split("samtools view -F 2316 -q {} {} {}".format(minMQ, fin_bam, region_value)))  # 2308  multiple threads -@ 4  is need ?
    ReadTree = list() 
    pre_start = 0; pre_end = 0; pre__seq = ""
    for l in p2.stdout:
        l = l.strip().split()
        if l[0][0] == "@":
            continue
        CIGAR = l[5]
        refStart = int(l[3]) - 1
        refEnd = refStart
        advance = 0
        for m in str(CIGAR):
            if m.isdigit():
                advance = advance * 10 + int(m)
                continue
            if m == "M" or m == "=" or m == "X" or m == "D":
                refEnd += advance
            advance = 0

        s = l[9].upper()
        
        # duplication?  else or not
        if refStart == pre_start and refEnd == pre_end and pre_seq == s:
            continue
        pre_start = refStart; pre_end = refEnd; pre_seq = s
        if useBaseQuality:
            ReadTree.append([refStart, refEnd, s, CIGAR, l[10]])
        else:
            ReadTree.append([refStart, refEnd, s, CIGAR])

    p2.stdout.close()
    p2.wait()
    
    ReadTree = sorted(ReadTree, key=lambda x: (x[0], x[1]))

    return ReadTree

def ksw_core(consensus, target, chrName, RNAME, reference_sequence, rS, shift, x_score, multiCan, max_merge_dis, Sig):
    canList = []
    for msa in consensus:
        if len(msa) == 0: continue
        #msa = bytes.decode(consensus[i])
        msa = bytes.decode(msa)
        alignment = ksw_aligner(msa, target, x_score)
        alignment[1] = getExactMatch(alignment[1], target, msa, 0)
        
        if Sig: print(alignment)
        L = makeCandidate(chrName, alignment, RNAME.split('_'), reference_sequence, rS, msa, shift, "LC", multiCan, max_merge_dis)
        canList.append(L)
    
    return canList

def overlap(s1, e1, s2, e2):
    S = max(s1, s2)
    E = min(e1, e2)

    if S < E: return True
    else: return False

def binarySearch(ReadTree, low, high, s, e, p):
    #s = p; e = p
    while low <= high:
        mid = (low + high) >> 1
        tmpl, tmpr = ReadTree[mid][0:2]
        if e < tmpl:
            high = mid - 1
        elif s > tmpr:
            low = mid + 1
        else:
            return mid

    return -1

def search_queries(ReadTree, sid, eid, S, E, p):
    queries = list()

    mid = binarySearch(ReadTree, sid, eid, S, E, p)
    bound = MAXR 
    A = sid if mid - bound < sid else mid - bound
    B = eid - 1 if mid + bound >= eid else mid + bound

    for i in range(A, B):
        t = ReadTree[i]
        if t[0] > E: break
        if t[1] < S: continue
        if overlap(t[0], t[1], S, E):
            queries.append(t)

    #for i in range(A, B):
    #    t = ReadTree[i]
    #    if t[0] > p: break
    #    if t[1] < p: continue
    #    if t[0] < p and t[1] > p:
    #        queries.append(t)

    return queries


def extractSeq(candidate, reference_sequence, queries, flanking, indL, refStart, useBaseQuality, useAllReads, minCNT, Sig):
    
    Start = candidate[0].pos - flanking
    End = candidate[-1].pos + flanking + indL

    query, qual, target, name = getLocalSeq(reference_sequence, refStart, queries, candidate, Start, End, useBaseQuality, useAllReads, flanking, Sig)

    if len(query) < minCNT: return None
    else: return [query, qual, target, name] 


def run_poa(candidate, reference_sequence, queries, refStart, flanking, useAllReads, indL, useBaseQuality, args, multiCan, IsID, Sig):
    # using abpoa for first round
    # flanking = 20; useAllReads = False
    # S = candidate[0].pos - flanking
    # E = candidate[-1].pos + flanking + indL
    res = extractSeq(candidate, reference_sequence, queries, flanking, indL, refStart, useBaseQuality, useAllReads, args.minCNT, Sig)
    if res is None: return [],[],[],[],""

    query, qual, target, name = res
    if Sig: print(name, target)

    res = POA(query, 1, args.ratio_identity_snp if IsID == "S" else args.ratio_identity_indel, is_MSA)
    if Sig: print("result of abPOA: reads cnt = ", len(query), res)
    canList = ksw_core(res[1], target, args.chrName, name, reference_sequence, refStart, args.shift, args.mismatch2 if multiCan else args.mismatch, multiCan, args.max_merge_dis, Sig)
    #if Sig: print(canList)

    return canList, res[0], query, qual, target

def run_ass(candidate, reference_sequence, queries, refStart, flanking, useAllReads, indL, useBaseQuality, multiCan, args, Sig):
    # using assemble for second round
    res = extractSeq(candidate, reference_sequence, queries, flanking, indL, refStart, useBaseQuality, useAllReads, args.minCNT, Sig)
    if res is None: return [],[],[],[], ""

    query, qual, target, name = res

    if Sig: print(name)
    rdl = [len(s) for s in query]
    minW = min(rdl) - 5
    maxW = max(rdl) + 5
    minW = minW if minW < 30 else 30
    maxW = maxW if maxW < 70 else 70
    res = mantaAss(query, minW, maxW)

    if Sig: print("result of manta: reads cnt = ", len(query), res)
    canList = ksw_core(res[1], target, args.chrName, name, reference_sequence, refStart, args.shift, args.mismatch2 if multiCan else args.mismatch, multiCan, args.max_merge_dis, Sig)

    return canList, res[0], query, qual, target


def main_ctrl(args):
    ################### check file exist ###########################
    CheckFileExist(args.fin_bam)
    CheckFileExist(args.fin_ref)
    CheckFileExist(args.fin_bam, sfx='.bai')
    CheckFileExist(args.fin_ref, sfx='.fai')
    # writing header
    vcfOut = open(args.fout_vcf, 'w')
    PrintVCFHeader(vcfOut, args.fin_ref, "HG002")
    ################################################################

    if args.target != -1:
        args.chrStart = args.target - 50000
        args.chrEnd = args.target + 50000
    
    ################### init values  ###############################
    useBaseQuality = args.useBaseQuality
    chrName = args.chrName
    refStart = args.chrStart
    refEnd = args.chrEnd
    is_range_given = args.chrStart is not None and args.chrEnd is not None
    if not is_range_given: 
        refStart = 1; refEnd = None

    P_ERROR = {"INS":args.perror_for_indel, "DEL":args.perror_for_indel, "S":args.perror_for_snp}
    ################################################################
    
    ###################### get ref, read, candidate ################
    reference_sequence = get_reference_sequence(args.fin_ref, args.chrName, refStart, refEnd)
    ReadTree = ExtractSeqAroundVariant(args.fin_bam, refStart, refEnd, args.chrName, args.minMQ, useBaseQuality)
    candidate_generator = candidate_from_file(args.fin_can, args.max_merge_dis, is_range_given, refStart, refEnd)
    ################################################################

    
    final_list = list()
    for candidate in candidate_generator:
        if not isinstance(candidate, list): #HC
            var = callHCVariant(candidate, P_ERROR)
            final_list.append(var)
        else:
            Sig = False
            for c in candidate:
                if c.pos == args.target: Sig = True

            multiCan = False if len(candidate) == 1 else True
            indL, ifLRP = getFlanking(candidate)
            IsID = "S" if indL == 0 else "ID"
            assty = "msa"
            
            ############# get total regional reads ##########################
            Start = candidate[0].pos - args.flanking
            End = candidate[-1].pos + args.flanking + indL
            Start = refStart if Start < refStart else Start
            End = refEnd if End > refEnd else End

            sid = 0; eid = len(ReadTree) - 1
            p = (End + Start) / 2
            queries = search_queries(ReadTree, sid, eid, Start, End, p)
            if len(queries) < args.minCNT: continue
            #################################################################

            if indL > 20: checksoft_and_realign(reference_sequence, refStart, queries, Start, End, indL, Sig) 
            
            ################ using abpoa for first round ####################
            flanking = 50; useAllReads = False
            canList, cluster_n, query, qual, target = run_poa(candidate, reference_sequence, queries, refStart, flanking, useAllReads, indL, useBaseQuality, args, multiCan, IsID, Sig)
            #################################################################

            ################ using assembly for second round ################
            if sum([len(s) for s in canList]) == 0 or args.poaORass == 2:
                flanking = args.flanking; useAllReads = True 
                canList, cluster_n, query, qual, target = run_ass(candidate, reference_sequence, queries, refStart, flanking, useAllReads, indL, useBaseQuality, multiCan, args, Sig)
                assty = "ass"

                if sum([len(s) for s in canList]) == 0: continue
                #if Sig: print(canList)
            #################################################################

            canList = mergeCoVar(canList, cluster_n, len(query), target, useBaseQuality, assty)

            if Sig:
                for c in canList:
                    print(c.pos, c.type, c.REF, c.ALT, c.supportA, c.supportT, c.con1, c.con2, c.read_T, c.assty)
            
            for c in canList:
                if c.type == "S":
                    most_likely_genotype([c], query, qual, args.mismatch2 if multiCan else args.mismatch, Sig)
                else:
                    callLCVariant([c], P_ERROR, -1)

            final_list.extend(canList)

    Generate_Output(vcfOut, final_list)
    vcfOut.close()

def run():
    parser = argparse.ArgumentParser(description="MSA examples using abPOA")
    parser.add_argument("--fin_bam", type=str, default="input_local.bam",
            help="Sorted bam file, default: %(default)s")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input")
    parser.add_argument("--fout_vcf", type=str, default="detect.vcf",
            help="vcf path for detected variants, default: %(default)s")
    parser.add_argument("--fin_can", type=str, default="PIPE",
            help="Candidate generated by ExtractVariantCandidate.py, default: %(default)s")
    parser.add_argument("--fin_repeats", type=str, default=None,
            help="Tandem repeat intervals for indel assembly, default: %(default)s")
    parser.add_argument("--chrName", type=str, default=None,
            help="The name of reference to be processed, default: %(default)s")
    parser.add_argument("--chrStart", type=int, default=None,
            help="The 1-based starting positions of the reference to be processed")
    parser.add_argument("--chrEnd", type=int, default=None,
            help="The inclusive ending positions of the reference to be processed")
    parser.add_argument("--perror_for_snp", type=float, default=0.1,
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")
    parser.add_argument("--perror_for_indel", type=float, default=0.1,
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")

    parser.add_argument("--ratio_identity_snp", type=float, default=0.2,
            help="min frequency of each consensus for heterozygote SNP, default:%(default)f")
    parser.add_argument("--ratio_identity_indel", type=float, default=0.2,
            help="min frequency of each consensus for heterozygote indel, default:%(default)f")
    parser.add_argument("--minMQ", type=int, default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")
    parser.add_argument("--minCNT", type=int, default=3,
            help="Minimum read counts required to call a variant, default:%(default)d")
    parser.add_argument("--max_merge_dis", type=int, default=5,
            help="Max distance to merge two variant candidates, default:%(default)d")
    parser.add_argument("--shift", type=int, default=5,
            help="The distance bewteen the detect variant and the activate region, default:%(default)d")
    parser.add_argument("--flanking", type=int, default=50,
            help="Flanking base pairs around variant site, default: %(default)d")
    parser.add_argument("--useBaseQuality", action='store_true', default=False,
            help="Use base quality to call variant, which will realignmet local reads to haplotypes and cost more time with little improvment of performance, default: True")
    parser.add_argument("--useAllReads", action='store_true', default=False,
            help="Use all local reads (including spanning reads and part-overlapped reads) to perform MSA, default: True")
    parser.add_argument("--target", type=int, default=-1,
            help=". default: %(default)d")
    parser.add_argument("--poaORass", type=int, default=1,
            help="POA or local assembly. default: %(default)d")
    parser.add_argument('-x', "--mismatch", type=int, default=4,
            help="Mismatch penalty. default: %(default)d")
    parser.add_argument('-d', "--mismatch2", type=int, default=6,
            help="Mismatch penalty. default: %(default)d")



    args = parser.parse_args()    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    t_s = time()
    main_ctrl(args)
    t_e = time()

    #logging.info("[localMSA] Finish region (%s:%d-%d) in %.2f seconds" %(args.chrName, args.chrStart, args.chrEnd, t_e - t_s))

if __name__ == "__main__":
    run()
