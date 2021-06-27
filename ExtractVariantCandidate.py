#!/usr/bin/env python
# coding=utf-8
import os
import sys
import argparse
import logging
import param
import shlex
from subprocess import PIPE
from collections import Counter
import gc
from time import time
from intervaltree import IntervalTree
from collections import defaultdict
from utils import setup_logging
from check import CheckFileExist, subprocess_popen

majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})
is_pypy = '__pypy__' in sys.builtin_module_names

SNP2S = {"A":"S", "a":"S", "C":"S", "c":"S", "G":"S", "g":"S", "T":"S", "t":"S", "I":"I", "D":"D", "N":"S", "n":"S"}
def divide_candidate(repeatMask, item, minCov_for_snp, minCov_for_indel, min_c_snp_HC, min_c_indel_HC, min_r_snp_HC, min_r_indel_HC, A, Ar, B, Br):
    chrN, pos, readC, refB, refC, var, alleleC = item[0:7]
    typ = var[0]
    if readC < 1:
        return "Filter"

    typ = SNP2S[typ] 
    
    if typ == "S":
        if alleleC < minCov_for_snp: return "Filter"
        if len(item) > 7: return "LC"  # multiallel   can alse be HC variant
        #else: return "HC"
        
        if repeatMask is not None and repeatMask.overlaps(pos):
            return "LRP"

        ratio = float(alleleC) / readC
        if readC >= min_c_snp_HC: #12
            return "HC" if ratio >= min_r_snp_HC else "LC"   #0.4, 0.3(best)
        elif readC >= A : # 5
            return "HC" if ratio >= Ar else "LC"
        else:
            return "LC"
    else:
        if len(item) > 7:
            var2, alleleC2 = item[7:9]
            if var2[0] != typ: return "LC"
            alleleC +=  alleleC2
        if alleleC == 1 or (alleleC < minCov_for_indel and float(alleleC)/readC < 0.75): return "Filter"


        if repeatMask is not None and repeatMask.overlaps(pos):
            return "LRP"

        ratio = float(alleleC) / readC
        if readC >= min_c_snp_HC: 
            return "HC" if ratio >= min_r_indel_HC else "LC"  #0.5
        elif readC >= B:
            return "HC" if ratio >= Br else "LC"
        else:
            return "LC"


class CandidateStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()

def overlap(s1, e1, s2, e2):
    S = max(s1, s2)
    E = min(e1, e2)

    if S < E: return True
    else: return False

def get_repeat_region(repeat_file_path, ctgName, ctgStart, ctgEnd):
    repeatMask = IntervalTree()
    with open(repeat_file_path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] != ctgName: continue
            s = int(line[1]); e = int(line[2])
            if overlap(s, e, ctgStart, ctgEnd):
                repeatMask.addi(s, e)

    return repeatMask

def get_reference_sequence(fasta_file_path, ctgName, ctgStart, ctgEnd):
    refernce_sequences = []
    if ctgStart != None and ctgEnd != None:
        region_value = "{}:{}-{}".format(ctgName, ctgStart, ctgEnd)
    else:
        region_value = str(ctgName)
    samtools_faidx_process = subprocess_popen(
        shlex.split("samtools faidx {} {}".format(fasta_file_path, region_value))
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


def MakeCandidates(args):
    CheckFileExist(args.fin_bam + ".bai")
    CheckFileExist(args.fin_ref + ".fai")
    chrName = args.chrName
    chrStart = args.chrStart
    chrEnd = args.chrEnd
    is_range_given = chrStart is not None and chrEnd is not None
    candidate_output_path = args.fout_can
    is_using_stdout_for_output_candidate = candidate_output_path == "PIPE"

    repeat = args.fin_repeat is not None

    if is_using_stdout_for_output_candidate:
        can_fp = CandidateStdout(sys.stdout)
    else:
        can_fpo = open(candidate_output_path, "wb")
        can_fp = subprocess_popen(shlex.split("gzip -c"), stdin=PIPE, stdout=can_fpo)
    
    if is_range_given:
        refStart = chrStart; refEnd = chrEnd
        refStart = 1 if refStart < 1 else refStart
    else:
        refStart = 1; refEnd = None
    #refStart -= param.expandReferenceRegion
    #refStart = 1 if refStart < 1 else refStart
    #refEnd += param.expandReferenceRegion

    RefSeq = get_reference_sequence(args.fin_ref, chrName, refStart, refEnd)
    if repeat:
        repeatMask = get_repeat_region(args.fin_repeat, chrName, refStart, refEnd)
    else:
        repeatMask = None

    Result_list = list() 
    pileup = defaultdict(lambda: {"A":0,"C":0,"D":0,"G":0,"I":0,"N":0,"T":0})
    #pileup_l = defaultdict(lambda: {"A":[], "C":[], "D":[], "G":[], "I":[], "N":[], "T":[]})
    pileup_l = defaultdict(lambda: {"A":0, "C":0, "D":[], "G":0, "I":[], "N":0, "T":0})
    
    if is_range_given:
        p2 = subprocess_popen(shlex.split("samtools view -F 2316 -q %d %s %s:%d-%d" % (args.minMQ, args.fin_bam, chrName, chrStart, chrEnd) ))   #2308   #  multiple threads -@ 4  is need ?
    else:
        p2 = subprocess_popen(shlex.split("samtools view -F 2316 -q %d %s %s" % (args.minMQ, args.fin_bam, chrName) ))

    A = args.min_c_snp_HC // 2
    Ar = args.min_r_snp_HC * 2 
    B = args.min_c_indel_HC // 2
    Br = args.min_r_indel_HC * 2
    Ar = Ar if Ar < 0.8 else 0.8
    Br = Br if Br < 0.8 else 0.8
    
    while True:
        row = p2.stdout.readline()
        is_finished = row == '' and p2.poll() is not None

        if row:
            l = row.strip().split()

            if l[0][0] == "@":
                continue

            QNAME = l[0]
            RNAME = l[2]

            if RNAME != chrName:
                continue
            
            FLAG = int(l[1])
            POS = int(l[3]) - 1 # switch from 1-base to 0-base to match sequence index
            MQ = int(l[4])
            CIGAR = l[5]
            SEQ = l[9].upper()
            QUAL = l[10]
            refPos = POS
            queryPos = 0

            advance = 0
            for c in str(CIGAR):
                if c.isdigit():
                    advance = advance * 10 + int(c)
                    continue

                if c == "S":
                    queryPos += advance
                elif c == "M" or c == "X" or c == "=":
                    for _ in range(advance):
                        tbase = SEQ[queryPos]
                        pileup[refPos][tbase] += 1
                        pileup_l[refPos][tbase] += ord(QUAL[queryPos]) - 33
                        refPos += 1
                        queryPos += 1
                elif c == "I":
                    if advance < 53:
                        pileup[refPos-1]["I"] += 1
                        pileup_l[refPos-1]["I"].append(SEQ[queryPos:queryPos+advance])
                    queryPos += advance
                elif c == "D":
                    if advance < 53:
                        pileup[refPos-1]["D"] += 1
                        pileup_l[refPos-1]["D"].append(RefSeq[refPos + 1 - refStart: refPos + 1 + advance - refStart])
                    refPos += advance
                advance = 0

        positions = [x for x in pileup.keys() if x < POS] if not is_finished else list(pileup.keys())
        positions.sort()
        for base_pos in positions:
            if is_range_given and not (chrStart <= base_pos + 1 <= chrEnd): continue
            flag = pileup.get(base_pos)
            if flag is None: continue


            baseCount = list(pileup[base_pos].items())
            baseCount_l = pileup_l[base_pos]
            refBase = RefSeq[base_pos + 1 - refStart]

            out = None
            out = CheckCandaidate(chrName, base_pos, baseCount, baseCount_l, refBase, args.minRatio_snp, args.minRatio_indel, args.variantType, args.min_ava_BQ)
            if out != None:
                hc_or_lc = divide_candidate(repeatMask, out, args.minCov_for_snp, args.minCov_for_indel, args.min_c_snp_HC, args.min_c_indel_HC, args.min_r_snp_HC, args.min_r_indel_HC, A, Ar, B, Br)
                if hc_or_lc != "Filter":
                    outp = " ".join([str(x) for x in out]) + " " + hc_or_lc + "\n"
                    can_fp.stdin.write(outp)
        for base_pos in positions:
            del pileup[base_pos]
        if is_finished:
            break

    p2.stdout.close()
    p2.wait()

    if not is_using_stdout_for_output_candidate:
        can_fp.stdin.close()
        can_fp.wait()
        can_fpo.close()

def CheckCandaidate(chrName, pos, baseCount, baseCount_l, refBase, minRatio_snp, minRatio_indel, variantType, baseQuality):
    if refBase == "N":
        return None
    totalCount = 0; readCount = 0; indelCount = 0; refCount = 0
    for x in baseCount:
        totalCount += x[1]
        if x[0] not in ["I", "D"]: 
            readCount += x[1]
            if x[0] == refBase:
                refCount += x[1]
        else:
            indelCount += x[1]

    denominator = totalCount
    if denominator == 0:
        denominator = 1
    baseCount.sort(key = lambda x:-x[1])  # sort baseCount descendingly
    
    #print(pos+1, totalCount, refBase, readCount, refCount, indelCount, baseCount, baseCount_l)

    p0 = float(baseCount[0][1]) / denominator
    
    sig = False
    if baseCount[0][0] == refBase and p0 > 1.0 - min(minRatio_snp, minRatio_indel):
        return None
    outC = [chrName, pos+1, readCount, refBase, refCount]
    for baseC in baseCount[:3]:
        p1 = float(baseC[1]) / denominator
        typ = baseC[0]
        if typ in ["I", "D"]: 
            threshold = minRatio_indel
        else:
            threshold = minRatio_snp
        if p1 >= threshold and baseC[0] != refBase:
            if typ in ["I", "D"]:  # indel candidate
                if variantType in ["all", "indel"]:
                    mostAllele = Counter(baseCount_l[typ]).most_common(2) 
                    for t in mostAllele:
                        outC.extend([typ+t[0],t[1]])
                    sig = True
            else: # SNP candidate
                if variantType in ["all", "snp"]:
                    ave_qual = baseCount_l[typ] / baseC[1]
                    if ave_qual > baseQuality:  # 20 best
                        outC.extend([typ, baseC[1]])
                        sig = True

    if sig: return outC

    return None

def main():
    parser = argparse.ArgumentParser(description="Generate variant candidates using alignments")
    parser.add_argument("--fin_bam", type=str, default="input_extract.bam",
            help="Sorted bam file, default: %(default)s")
    parser.add_argument("--fin_ref", type=str, default="ref.fa",
            help="Reference fasta input, default: %(default)s")
    parser.add_argument("--fout_can", type=str, default="PIPE",
            help="Variant candidate output prefix, default: %(default)s")
    parser.add_argument("--fin_repeat", type=str, default=None,
            help="Tandem repeat intervals for indel assembly, default: %(default)s")
    parser.add_argument("--minMQ", type=int, default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")
    parser.add_argument("--min_ava_BQ", type=int, default=20,
            help="Minimum average base quality required to consider a base for calling, default:%(default)d")
    parser.add_argument("--minCov_for_snp", type=int, default=2,
            help="Minimum read counts required to call a snp, default:%(default)d")
    parser.add_argument("--minCov_for_indel", type=int, default=4,
            help="Minimum read counts required to call a indel, default:%(default)d")
    parser.add_argument("--minRatio_snp", type=float, default=0.125,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for SNP, default:%(default)f")
    parser.add_argument("--minRatio_indel", type=float, default=0.1,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for INDEL, default:%(default)f")
    parser.add_argument("--min_r_snp_HC", type=float, default=0.3,
            help="Min  supporting reads count ratio to detect a high-confidence SNP variant candidate, default:%(default)f")
    parser.add_argument("--min_r_indel_HC", type=float, default=0.5,
            help="Min  supporting reads count ratio to detect a high-confidence indel variant candidate, default:%(default)f")
    parser.add_argument("--min_c_snp_HC", type=int, default=12,
            help="Min  supporting reads count to detect a high-confidence SNP variant candidate, default:%(default)d")
    parser.add_argument("--min_c_indel_HC", type=int, default=10,
            help="Min  supporting reads count to detect a high-confidence indel variant candidate, default:%(default)d")
    parser.add_argument("--chrName", type=str, default="22",
            help="The name of reference to be processed, default:%(default)s")
    parser.add_argument("--chrStart", type=int, default=None,
            help="The 1-based starting positions of the reference to be processed")
    parser.add_argument("--chrEnd", type=int, default=None,
            help="The inclusive ending positions of the reference to be processed")
    parser.add_argument("--variantType", type=str, default="all",
            help="Extract candidates of SNP, indel or all of small variant,[snp, indel, all], default: %(default)s")
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()

    t_s = time()
    MakeCandidates(args)
    t_e = time()
    #logging.info("[ExtractVariantCandidate] Finish region (%s:%d-%d) in %.2f seconds" %(args.chrName, args.chrStart, args.chrEnd, t_e - t_s))

if __name__ == "__main__":
    main()
