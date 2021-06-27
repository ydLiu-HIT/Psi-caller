#!/usr/bin/env python
# coding=utf-8
import sys
import shlex
import logging
import argparse
import signal
from time import time
from utils import setup_logging
from os.path import dirname

from command_options import *
from check import CheckCmdExist, CheckFileExist, subprocess_popen

class InstancesClass(object):
    def __init__(self):
        self.extract_variant_candidate = None
        self.local_MSA = None

    def poll(self):
        self.extract_variant_candidate.poll()
        self.local_MSA.poll()


c = InstancesClass()


def check_return_code(signum, frame):
    c.poll()
    #print >> sys.stderr, c.extract_variant_candidate.returncode, c.create_tensor.returncode, c.call_variant.returncode
    if c.extract_variant_candidate.returncode != None and c.extract_variant_candidate.returncode != 0:
        c.local_MSA.kill()
        sys.exit("ExtractVariantCandidates.py exited with exceptions. Exiting...")

    if c.local_MSA.returncode != None and c.local_MSA.returncode != 0:
        c.extract_variant_candidate.kill()
        sys.exit("localMSA.py exited with exceptions. Exiting...")

    if (
        c.extract_variant_candidate.returncode == None or
        c.local_MSA.returncode == None 
    ):
        signal.alarm(5)


def run(args):
    basedir = dirname(__file__)
    EVCBin = basedir + "ExtractVariantCandidate.py"
    LMSABin = basedir + "localMSA.py"
    pypyBin = CheckCmdExist(args.pypy)

    fin_bam = args.fin_bam
    fin_ref = args.fin_ref
    fin_repeat = args.fin_repeat

    chrName = args.chrName
    variantType = args.variantType
    minMQ = args.minMQ
    minBQ = args.min_ava_BQ
    minRatio_snp = args.minRatio_snp
    minRatio_indel = args.minRatio_indel
    minCNT = args.minCNT
    shift = args.shift
    flanking = args.flanking
    perror_for_snp = args.perror_for_snp
    perror_for_indel = args.perror_for_indel
    ratio_identity_snp = args.ratio_identity_snp
    ratio_identity_indel = args.ratio_identity_indel
    
    minCov_for_snp = args.minCov_for_snp
    minCov_for_indel = args.minCov_for_indel
    min_r_snp_HC = args.min_r_snp_HC
    min_r_indel_HC = args.min_r_indel_HC
    min_c_snp_HC = args.min_c_snp_HC
    min_c_indel_HC = args.min_c_indel_HC
    max_merge_dis = args.max_merge_dis

    region_chunk_size = args.chunkWidth
    output_prefix = args.fout_vcf

    chrStart = None
    chrEnd = None
    if args.chrStart is not None and args.chrEnd is not None and int(args.chrStart) <= int(args.chrEnd):
        chrStart = CommandOption('chrStart', args.chrStart)
        chrEnd = CommandOption('chrEnd', args.chrEnd)

    #Exec
    extract_variant_command_options = [
        pypyBin,
        EVCBin,
        CommandOption('fin_bam', fin_bam),
        CommandOption('fin_ref', fin_ref),
        CommandOption('fin_repeat', fin_repeat),
        CommandOption('minMQ', minMQ),
        CommandOption('min_ava_BQ', minBQ),
        CommandOption('minRatio_snp', minRatio_snp),
        CommandOption('minRatio_indel', minRatio_indel),
        CommandOption("variantType", variantType),
        CommandOption('minCov_for_snp', minCov_for_snp),
        CommandOption('minCov_for_indel', minCov_for_indel),
        CommandOption('min_r_snp_HC', min_r_snp_HC),
        CommandOption('min_r_indel_HC', min_r_indel_HC),
        CommandOption('min_c_snp_HC', min_c_snp_HC),
        CommandOption('min_c_indel_HC', min_c_indel_HC),
        CommandOption('chrName', chrName),
        chrStart,
        chrEnd
    ]

    local_MSA_command_options = [
        pypyBin,
        LMSABin,
        CommandOption('fin_bam', fin_bam),
        CommandOption('fin_ref', fin_ref),
        CommandOption('minMQ', minMQ),
        CommandOption('perror_for_snp', perror_for_snp),
        CommandOption('perror_for_indel', perror_for_indel),
        CommandOption('ratio_identity_snp', ratio_identity_snp),
        CommandOption('ratio_identity_indel', ratio_identity_indel),
        CommandOption('shift', shift),
        CommandOption('flanking', flanking),
        CommandOption("max_merge_dis", max_merge_dis),
        CommandOption('chrName', chrName),
        chrStart,
        chrEnd,
        CommandOption("fout_vcf", output_prefix),
    ]

    if args.useBaseQuality:
        local_MSA_command_options.append(CommandOptionWithNoValue('useBaseQuality'))

    if args.useAllReads:
        local_MSA_command_options.append(CommandOptionWithNoValue('useAllReads'))


    try:
        c.extract_variant_candidate = subprocess_popen(
            shlex.split(command_string_from(extract_variant_command_options))
        )

        c.local_MSA = subprocess_popen(
            shlex.split(command_string_from(local_MSA_command_options)),
            stdin=c.extract_variant_candidate.stdout
        )
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit("Failed to start required processes. Exiting...")

    signal.signal(signal.SIGALRM, check_return_code)
    signal.alarm(2)

    try:
        c.local_MSA.wait()
        c.extract_variant_candidate.stdout.close()
        c.extract_variant_candidate.wait()
    except KeyboardInterrupt as e:
        print("KeyboardInterrupt received when waiting at calling.py, terminating all scripts.")
        try:
            c.local_MSA.terminate()
            c.extract_variant_candidate.terminate()
        except Exception as e:
            print(e)

        raise KeyboardInterrupt
    except Exception as e:
        print("Exception received when waiting at calling.py, terminating all scripts.")
        print(e)
        try:
            c.local_MSA.terminate()
            c.extract_variant_candidate.terminate()
        except Exception as e:
            print(e)

        raise e


def main():
    parser = argparse.ArgumentParser(prog="calling",
            description="Call variant using sorted bam",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    
    #*******Parameters of ExtractVariantCandidate.py****
    parser.add_argument("--fin_bam", 
            type=str,
            default="input_calling.bam",
            help="Sorted bam file, [BAM].")

    parser.add_argument("--fin_ref",
            type=str,
            default="ref.fa",
            help="Reference fasta input, [FASTA]")

    parser.add_argument("--fout_vcf", 
            type=str, 
            default="detect.vcf",
            help="vcf path for detected variants, default: %(default)s")

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

    parser.add_argument("--minRatio_snp", 
            type=float,
            default=0.125,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for SNP, default:%(default)f")

    parser.add_argument("--minRatio_indel", 
            type=float,
            default=0.1,   ## previouse use 0.125 as default
            help="Minimum variant supported read count ratio for INDEL, default:%(default)f")

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

    parser.add_argument("--chrName", 
            type=str,
            default=None,   ## previouse use 0.125 as default
            help="Contig name of sequence to be processed")
    
    parser.add_argument("--chrStart", 
            type=int,
            default=None,   ## previouse use 0.125 as default
            help="The 1_based starting position of sequence to be processed")

    parser.add_argument("--chrEnd", 
            type=int,
            default=None,   ## previouse use 0.125 as default
            help="The 1_based inclusive ending position of sequence to be processed")

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
    
    setup_logging()
    t_s = time()
    run(args)
    t_e = time()
    logging.info("[calling.py] Finish region (%s:%d-%d) in %.2f seconds" %(args.chrName, args.chrStart, args.chrEnd, t_e - t_s))


if __name__ == "__main__":
    main()
