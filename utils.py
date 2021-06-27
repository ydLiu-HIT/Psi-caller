#!/usr/bin/env python
# coding=utf-8

import logging
from check import CheckFileExist

#majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})

base_to_ACGT = dict(zip(
    "ACGTURYSWKMBDHVN",
    ("A", "C", "G", "T", "T", "A", "C", "C", "A", "G", "A", "C", "A", "A", "A", "A")
))

def PrintVCFHeader(file, fin_ref, sample):
    contigINFO = list() 
    fai_fn = CheckFileExist(fin_ref + ".fai")

    with open(fai_fn, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] in majorContigs:
                contigINFO.append([line[0], int(line[1])])

    file.write('##fileformat=VCFv4.1\n')

    for i in contigINFO:
        file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

    file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    file.write('##FILTER=<ID=FILTER,Description="Genotyping model thinks this site is reference or lower than a threshold score.">\n')
    file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    file.write('##INFO=<ID=LENGUESS,Number=1,Type=Integer,Description="Best guess of the indel length">\n')
    file.write('##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Consensus id of MSA">\n')
    file.write('##INFO=<ID=MQ,Number=60,Type=Integer,Description="Mapping quality of MSA">\n')
    file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
    file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read Depth for each allele">\n')
    file.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1)">\n')


    file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (sample))

def Generate_Output(file, Result):
    idx = 0
    prePos = 0
    preTyp = "INS"
    preAlt = ""
    for item in Result:
        GQ = item.GQ
        item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.con1, item.con2, item.read_T]
        if item[1] == prePos and item[2] == preTyp and item[4] == preAlt:
            prePos = item[1]
            preTyp = item[2]
            preAlt = item[4]
            continue
        info_list = "SVTYPE={SVTYPE};LENGUESS={LENGUESS};CONSENSUS={CON}".format(
            SVTYPE = item[2],
            LENGUESS = item[8], ##abs(len(item[3]) - len(item[4])),
            CON = item[5])
        #Cn = 0
        #for i in str(item[10]).split(','):
        #    Cn += int(i)
        #AD = str(item[10]) + "," + str(item[11]-Cn)
        AD = str(item[12])  + "," + str(item[13])
        file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}\n".format(
            CHR = item[0],
            POS = item[1],
            ID = ".",   #ID = "MLSV.%d"%(idx),
            REF = item[3],
            ALT = item[4],
            QUAL = item[9], #item[9],
            FILTER = "FILTER" if item[7] == "0/0" else "PASS",
            INFO = info_list,
            FORMAT = "GT:GQ:DP:AD",
            SAMPLE = item[7] + ":" + str(GQ) + ":" + str(item[14]) + ":" + AD))
        prePos = item[1]
        preTyp = item[2]
        preAlt = item[4]
        idx += 1

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)


