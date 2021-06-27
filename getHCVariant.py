#!/usr/bin/env python
# coding=utf-8
import sys, os
from time import time
from collections import Counter
import logging
import argparse
from multiprocessing import Pool
from Canclass import *
from math_func import *
from utils import Generate_Output

majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def filter_majorContig(fp):
    fin = open(fp, 'r')
    Lines = list()
    for line in fin:
        if line.split()[0] not in majorContigs:
            continue
        Lines.append(line.strip())
    fin.close()

    return Lines

def callHCVariant(candidate, P_ERROR):
    line = candidate.strip().split(' ')
    chrN, pos, readC, refB, refC, sv, alleleC = line[0:7]
    svT = sv[0]
    if svT == "I" or svT == "D":
        svB = sv[1:]
        Rs = refB+svB if svT=="D" else refB
        As = refB if svT=="D" else refB+svB
        if len(line) > 8:
            alleleC = int(alleleC) + int(line[8])
        item = [chrN,pos,"INS" if svT=="I" else "DEL",Rs,As,line[-1],60,len(svB),int(alleleC),int(readC), 0, 0]
        refC, readC = rescale_read_counts(item[9]-item[8], item[9])
    else:
        item = [chrN,pos,"S",refB,svT,"HC",60,0,int(alleleC),int(readC), 0, 0]
        refC, readC = rescale_read_counts(int(refC), item[9])
    var = variantsCan(item, 0, 0, 0, "HC")
    predictions = cal_GL(refC, readC, P_ERROR[item[2]])
    var.GT, var.GQ, var.qual = computer_quals(predictions, 2)

    return var

def callLCVariant(candidate, P_ERROR, target):
    multip = True if len(candidate) > 1 else False
    for var in candidate:
        Sig = False
        if var.pos == target: Sig = True
        svT = var.type.split(',')
        if len(svT) == 1:
            if multip:
                refC, readC = rescale_read_counts(var.con2, var.con2 + var.con1)
                #refC, readC = rescale_read_counts(var.supportT - var.supportA, var.supportT)
            else:
                if var.assty == "msa":
                    refC, readC = rescale_read_counts(var.con2, var.con1 + var.con2)
                elif var.assty == "ass":
                    #refC, readC = rescale_read_counts(var.read_T - var.con1, var.read_T)
                    refC, readC = rescale_read_counts(var.con2, var.con1 + var.con2)
                #refC, readC = rescale_read_counts(var.supportT - var.supportA, var.supportT)
            
            predictions = cal_GL(refC, readC, P_ERROR[svT[0]])
            var.GT, var.GQ, var.qual = computer_quals(predictions, 2)
        else:
            which = 1 if var.con1 < var.con2 else 0
            alleleC = var.con1 if which == 0 else var.con2
            refC, readC = rescale_read_counts(var.read_T - alleleC, var.read_T)

            predictions = cal_GL(refC, readC, P_ERROR[svT[which]])
            var.GT, var.GQ, var.qual = computer_quals(predictions, 2)

            if var.GT == "0/1":
                var.GT = "1/2"
                var.supportA = str(var.con1)+","+str(var.con2)
            elif var.GT == "1/1":
                var.type = svT[which]
                var.ALT = var.ALT.split(',')[which]
                var.varlen = var.varlen.split(',')[which]


#def get_hc_variantList(fp, threads, perror_for_indel, perror_for_snp):
#    Lines = filter_majorContig(fp)
#    hcCount = len(Lines)
#    analysis_pools = Pool(processes=int(threads))
#    chunks = hcCount//int(threads)
#    HCVariant = list(); semiResult = list()
#    start = 0; end = chunks
#    print(threads)
#    for m in range(int(threads)):
#        param = [(Lines[start:end], perror_for_indel, perror_for_snp)]
#        semiResult.append(analysis_pools.map_async(run_core, param))
#        #param = (Lines[start:end], perror_for_indel, perror_for_snp)
#        #HCVariant.extend(run_core(param))
#        start = end; end += chunks
#    #last
#    if start < hcCount:
#        param = [(Lines[start:hcCount], perror_for_indel, perror_for_snp)]
#        semiResult.append(analysis_pools.map_async(run_core, param))
#        #param = (Lines[start:hcCount], perror_for_indel, perror_for_snp)
#        #HCVariant.extend(run_core(param))
#
#    analysis_pools.close()
#    analysis_pools.join()
#    for res in semiResult:
#        try: HCVariant.extend(res.get()[0])
#        except: pass
#
#    del Lines; del semiResult
#
#    logging.info("Call HC variants finished")
#    return HCVariant
#
#def run_core(args):
#    return getVariants(*args)

#def getVariants(Vlist, perror_for_indel, perror_for_snp):
#    P_ERROR = {"INS":perror_for_indel, "DEL":perror_for_indel, "S":perror_for_snp}
#    CandidateList = []
#    for line in Vlist:
#        line = line.strip().split(' ')
#        chrN, pos, readC, refB, refC, sv, alleleC = line[0:7]
#        svT = sv[0]
#        if svT == "I" or svT == "D":
#            svB = sv[1:]
#            Rs = refB+svB if svT=="D" else refB
#            As = refB if svT=="D" else refB+svB
#            refC = int(readC) - int(alleleC)
#            item = [chrN,pos,"INS" if svT=="I" else "DEL",Rs,As,"HC",60,len(svB),refC,int(readC),0,0]
#        else:
#            item = [chrN,pos,"S",refB,svT,"HC",60,0,int(refC),int(readC),0,0]
#        var = variantsCan(item)
#        refC, readC = rescale_read_counts(item[8], item[9])
#        log10_probs = calc_reference_confidence(refC, readC, P_ERROR[item[2]])
#        real_probs = toRealSpace(log10_probs)
#        predictions = round_gls(real_probs)
#        var.GQ, var.qual, gt = computer_quals(predictions, 2)
#        var.GT = str(gt[0]) + '/' + str(gt[1]) 
#        if var.GT != "0/0": var.qual += 5.0
#
#        CandidateList.append(var)
#            
#    return CandidateList
#
#def get_hc_variant(args):
#    P_ERROR = {"INS":args.perror_for_indel, "DEL":args.perror_for_indel, "S":args.perror_for_snp}
#    fin = open(args.fin_info, 'r')
#    vcfout = open("tmp.vcf", 'w')
#    CandidateList = []
#    for line in fin:
#        line = line.strip().split(' ')
#        chrN, pos, readC, refB, refC, sv, alleleC = line[0:7]
#        svT = sv[0]
#        if svT == "I" or svT == "D":
#            svB = sv[1:]
#            Rs = refB+svB if svT=="D" else refB
#            As = refB if svT=="D" else refB+svB
#            refC = int(readC) - int(alleleC)
#            item = [chrN,pos,"INS" if svT=="I" else "DEL",Rs,As,"HC",60,len(svB),refC,int(readC),0,0]
#        else:
#            item = [chrN,pos,"S",refB,svT,"HC",60,0,int(refC),int(readC),0,0]
#
#        var = variantsCan(item)
#        refC, readC = rescale_read_counts(item[8], item[9])
#        log10_probs = calc_reference_confidence(refC, readC, P_ERROR[item[2]])
#        real_probs = toRealSpace(log10_probs)
#        predictions = round_gls(real_probs)
#        var.GQ, var.qual, gt = computer_quals(predictions, 2)
#        var.GT = str(gt[0]) + '/' + str(gt[1]) 
#        if var.GT != "0/0": var.qual += 5.0
#        CandidateList.append(var)
#        
#        #print predictions
#        #print var.GQ, var.qual, var.GT
#            
#    fin.close()
#    Generate_Output(vcfout, CandidateList)
#    vcfout.close()
#
#    cmd = "cat %s %s > %s" %("header", "tmp.vcf", args.fout_vcf)
#    os.system(cmd)


def run():
    parser = argparse.ArgumentParser(description="Get high confidence variant from CIAGR signals directly")
    parser.add_argument("--fin_info", type=str, default="hc.info",
            help="High confidence variant candidate generated by ExtractVariant.py, default: %(default)s")
    parser.add_argument("--fout_vcf", type=str, default="hc.vcf",
            help="High confidence variants output, default: %(default)s")
    parser.add_argument("--perror_for_snp", type=float, default=0.1,
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")
    parser.add_argument("--perror_for_indel", type=float, default=0.1,
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")
    parser.add_argument("--threads", type=int, default=1,
            help="Number of threads to use, default:%(default)d")
    args = parser.parse_args()    

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    get_hc_variant(args)
    #get_hc_variantList(args.fin_info, args.threads, args.perror_for_indel, args.perror_for_snp)


if __name__ =='__main__':
    t_s = time()
    run()
    logging.info("Finish the script in %f seconds", time() - t_s)
