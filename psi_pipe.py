#!/usr/bin/env python
# coding=utf-8
''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  psi-pipe.py
 * @author: Yadong Liu
 * @date: Jun 24th 2020
 * @version V1.0   
'''
import os
import sys
import argparse

from command_options import *
from check import CheckCmdExist, CheckFileExist

majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})

VERSION = '1.0'
PROG = "Psi-caller"

class maindp(object):
	'''
	Detailed descriptions of Psi-caller version and its parameters.

	'''

	USAGE="""\
	Fast and accurate SNP and indels detection for short Reads with %s.
	
	Current version: v%s
	Author: Yadong Liu
	Contact: ydliu@hit.edu.cn

	"""%(PROG, VERSION)

def run(args):
    basedir = os.path.dirname(__file__)

    callBin = basedir + "calling.py"
    pypyBin = CheckCmdExist(args.pypy)

    fin_bam = CheckFileExist(args.fin_bam)
    fin_bam_bai = CheckFileExist(fin_bam, sfx=".bai")
    fin_ref = CheckFileExist(args.fin_ref)
    fin_ref_fai = CheckFileExist(fin_ref, sfx=".fai")
    fin_repeat = args.fin_repeat
    if fin_repeat is not None:
        fin_repeat = CheckFileExist(args.fin_repeat)

    minMQ = args.minMQ
    minBQ = args.min_ava_BQ
    variantType = args.variantType
    minRatio_snp = args.minRatio_snp
    minRatio_indel = args.minRatio_indel
    minCNT = args.minCNT
    shift = args.shift
    flanking = args.flanking
    minCov_for_snp = args.minCov_for_snp
    minCov_for_indel = args.minCov_for_indel
    min_r_snp_HC = args.min_r_snp_HC
    min_r_indel_HC = args.min_r_indel_HC
    min_c_snp_HC = args.min_c_snp_HC
    min_c_indel_HC = args.min_c_indel_HC
    max_merge_dis = args.max_merge_dis

    perror_for_snp = args.perror_for_snp
    perror_for_indel = args.perror_for_indel
    ratio_identity_snp = args.ratio_identity_snp
    ratio_identity_indel = args.ratio_identity_indel
    region_chunk_size = args.chunkWidth

    if not os.path.exists(args.workDir):
        os.makedirs(args.workDir)
    output_prefix = args.workDir + '/' if args.workDir[-1] != '/' else args.workDir

    callVar_command_options = [
        ExecuteCommand('python', callBin),
        CommandOption('fin_bam', fin_bam),
        CommandOption('fin_ref', fin_ref),
        CommandOption('fin_repeat', fin_repeat),
        CommandOption('minMQ', minMQ),
        CommandOption('min_ava_BQ', minBQ),
        CommandOption('minCNT', minCNT),
        CommandOption('minRatio_snp', minRatio_snp),
        CommandOption('minRatio_indel', minRatio_indel),
        CommandOption('pypy', pypyBin),
        CommandOption("variantType", variantType),
        CommandOption('minCov_for_snp', minCov_for_snp),
        CommandOption('minCov_for_indel', minCov_for_indel),
        CommandOption('min_r_snp_HC', min_r_snp_HC),
        CommandOption('min_r_indel_HC', min_r_indel_HC),
        CommandOption('min_c_snp_HC', min_c_snp_HC),
        CommandOption('min_c_indel_HC', min_c_indel_HC),
        CommandOption('perror_for_snp', perror_for_snp),
        CommandOption('perror_for_indel', perror_for_indel),
        CommandOption('ratio_identity_snp', ratio_identity_snp),
        CommandOption('ratio_identity_indel', ratio_identity_indel),
        CommandOption("max_merge_dis", max_merge_dis),
        CommandOption('shift', shift),
        CommandOption('flanking', flanking),
    ]
    
    if args.useBaseQuality:
        callVar_command_options.append(CommandOptionWithNoValue('useBaseQuality'))

    if args.useAllReads:
        callVar_command_options.append(CommandOptionWithNoValue('useAllReads'))

    with open(fin_ref_fai, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")

            contig_name = columns[0]
            if str(contig_name) not in majorContigs:
                continue

            region_start, region_end = 1, 1
            contig_length = int(columns[1])
            while region_end < contig_length:
                region_start = region_end
                region_end = region_start + region_chunk_size
                if region_end > contig_length:
                    region_end = contig_length

                output_fn = "%s.%s_%d_%d.vcf" % (output_prefix+"var", contig_name, region_start, region_end)

                additional_command_options = [
                    CommandOption('chrName', contig_name),
                    CommandOption('chrStart', region_start),
                    CommandOption('chrEnd', region_end),
                    CommandOption('fout_vcf', output_fn)
                ]
                
                print(command_string_from(callVar_command_options) + " " + command_string_from(additional_command_options))

def main():
    parser = argparse.ArgumentParser(prog="Psi-caller",
            description=maindp.USAGE,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--version', '-v', 
		action = 'version', 
		version = '%(prog)s {version}'.format(version=VERSION))

    #*******Parameters of ExtractVariantCandidate.py****
    parser.add_argument("fin_bam", 
            type = str,
            help="Sorted bam file, [BAM].")

    parser.add_argument("fin_ref",
            type=str, 
            help="Reference fasta input, [FASTA]")

    parser.add_argument("workDir",
            type=str,
            help="Work-directory for distributed job")

    parser.add_argument("--fin_repeat", 
            type=str, 
            default=None,
            help="Tandem repeat intervals for indel assembly, default: %(default)s")

    parser.add_argument("--minMQ", 
            type=int, 
            default=10,
            help="Minimum Mapping Quality. Mapping quality lower than the setting will be filtered, default:%(default)d")

    parser.add_argument("--min_ava_BQ", 
            type=int, 
            default=20,
            help="Minimum average base quality required to consider a base for calling, default:%(default)d")

    parser.add_argument("--minCov_for_snp", 
            type=int, 
            default=2,
            help="Minimum read counts required to call a snp, default:%(default)d")

    parser.add_argument("--minCov_for_indel",
            type=int,
            default=4,
            help="Minimum read counts required to call a indel, default:%(default)d")

    parser.add_argument("--min_r_snp_HC", 
            type=float, 
            default=0.3,
            help="Min supporting reads count ratio to detect a high-confidence SNP variant candidate, default:%(default)f")

    parser.add_argument("--min_r_indel_HC", 
            type=float, 
            default=0.5,
            help="Min supporting reads count ratio to detect a high-confidence indel variant candidate, default:%(default)f")

    parser.add_argument("--min_c_snp_HC", 
            type=int, 
            default=12,
            help="Min supporting reads count to detect a high-confidence SNP variant candidate, default:%(default)d")

    parser.add_argument("--min_c_indel_HC", 
            type=int, 
            default=10,
            help="Min supporting reads count to detect a high-confidence indel variant candidate, default:%(default)d")

    parser.add_argument("--minRatio_snp", 
            type=float,
            default=0.125,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for SNP, default:%(default)f")

    parser.add_argument("--minRatio_indel", 
            type=float,
            default=0.1,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for INDEL, default:%(default)f")

    parser.add_argument("--variantType", 
            choices = ["all", "snp", "indel"], 
            type = str,
            default="all",
            help="Extract candidates of SNP, indel or all of small variant, default: %(default)s")

    parser.add_argument('-t', "--threads",
            type=int, 
            default=1,
            help="Number of threads to use, default:%(default)d")

    parser.add_argument("--chunkWidth", 
            type=int, 
            default=10000000,
            help="Reference length to detect candidate in one loop, default:%(default)d")

    parser.add_argument("--pypy", 
            type=str, 
            default="pypy3",
            help="Path to pypy, default:%(default)s")

    #*******Parameters of ExtractVariantCandidate.py****
    parser.add_argument("--perror_for_snp", 
            type=float, 
            default=0.1,
            help="P-error is the probability of observing a heterozygote SNP, default:%(default)f")

    parser.add_argument("--perror_for_indel", 
            type=float,
            default=0.1,
            help="P-error is the probability of observing a heterozygote indel, default:%(default)f")

    parser.add_argument("--ratio_identity_snp", 
            type=float, 
            default=0.2,
            help="min frequency of each consensus for heterozygote SNP, default:%(default)f")

    parser.add_argument("--ratio_identity_indel", 
            type=float, 
            default=0.2,
            help="min frequency of each consensus for heterozygote indel, default:%(default)f")

    parser.add_argument("--minCNT", 
            type=int, 
            default=3,
            help="Minimum read counts required to call a variant, default:%(default)d")

    parser.add_argument("--max_merge_dis", 
            type=int, 
            default=5,
            help="Max distance to merge two variant candidates, default:%(default)d")
    
    parser.add_argument("--shift", 
            type=int, 
            default=5,
            help="The distance bewteen the detect variant and the activate region, default:%(default)d")

    parser.add_argument("--flanking", 
            type=int, 
            default=50,
            help="Flanking base pairs around variant site which need to be extract for multiple sequence alignment, default: %(default)d")

    parser.add_argument("--useAllReads", 
            action='store_true', 
            default=False,
            help="Use all local reads (including spanning reads and part-overlapped reads) to perform MSA, default: True")

    parser.add_argument("--useBaseQuality",
            action='store_true', 
            default=False,
            help="Use base quality to call variant, which will realignmet local reads to haplotypes and cost more time with little improvment of performance, default: True")



    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    run(args)

if __name__ == "__main__":
    main()
