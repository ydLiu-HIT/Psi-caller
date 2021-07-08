# Psi-caller
a lightweight short read-based variant caller with high speed and accuracy



## Getting Start
    git clone https://github.com/ydLiu-HIT/Psi-caller.git
    cd Psi-caller/library/
    cd ksw2/; pypy3 setup.py install
    cd abpoa/; pypy3 setup.py install
    cd mantalib/; pypy3 setup.py install
    
    ## generate subtasks and variant calling
    python generate_task.py aln.bam reference.fa workspace/ > tasks.sh
    cat tasks.sh | parallel -j 6
    bash post_process.sh workspace output.vcf

## Introduction

With the rapid development of short read sequencing technologies, a number of population-scale resequencing studies have been carried out to study the associations between human genome variants and various phenotypes in recent years. Variant calling is one of the core bioinformatics tasks in such studies to comprehensively discover genomic variants in sequenced samples. Many efforts have been made to develop short read-based variant calling approaches, however, state-of-the-art tools are still computationally expensive, meanwhile, cutting-edge genomics studies also have higher requirements on the yields of variant calling. Herein, we propose Partial Order Alignment-based SNV and indel caller (Psi-caller), a lightweight variant calling algorithm that simultaneously achieves high performance and yield. Mainly, Psi-caller recognizes and divides the candidate variant site into three categories according to the complexity and location of the signatures, and employs various methods including binomial model, partial-order alignment and de-Bruijn graph-based local assembly to handle various categories of candidate sites to call and genotype SNVs/Indels, respectively. Benchmarks on simulated and real short-read sequencing datasets demonstrate that Psi-caller is times faster than state-of-the-art tools with higher or equal sensitivity and accuracy. It has the potential to well-handle large-scale datasets in cutting-edge genomics studies.

## Dependence

```
python3
pypy3
Cython
numpy
vcfcat
intervaltree
samtools  # samtools was used to extract read sequence in order to generate allele sequence.

### install dependence
conda install -c conda-forge pypy3.6
conda install -c bioconda vcflib
pypy3 -m ensurepip
pypy3 -m pip install Cython intervaltree numpy
```

## Synopsis

```
## extract variant candidate sites from pileups in the region (1:1-20000001)
pypy3 ExtractVariantCandidate.py --fin_bam aln.bam --fin_ref hs37d5.fa --chrName 1 --chrStart 1 --chrEnd 20000001 --fout_can workspace/var.1_1_20000001.can

## variant calling and genotyping from the candidates
pypy3 localMSA.py --fin_bam aln.bam --fin_ref hs37d5.fa --chrName 1 --chrStart 1 --chrEnd 20000001 --fin_can workspace/var.1_1_20000001.can --fout_vcf workspace/var.1_1_20000001.vcf

## generate subtasks and variant calling using multiple threads
python generate_task.py aln.bam reference.fa workspace/ > tasks.sh
cat tasks.sh | parallel -j 6

## generate separeted subtasks for Candidate recognition and variant calling using multiple threads
python separeted_task.py aln.bam reference.fa workspace/ --task_prefix subtask
cat subtask_extract.sh | parallel -j 6  ## Candidate recognition
cat subtask_call.sh | parallel -j 6  ## Variants calling 
```

## Commands and options

**ExtractVariantCandidate.py**

```
  -h, --help            show this help message and exit
  --fin_bam FIN_BAM     Sorted bam file, default: input_extract.bam
  --fin_ref FIN_REF     Reference fasta input, default: ref.fa
  --fout_can FOUT_CAN   Variant candidate output prefix, default: PIPE
  --fin_repeat FIN_REPEAT
                        Tandem repeat intervals for indel assembly, default:
                        None
  --minMQ MINMQ         Minimum Mapping Quality. Mapping quality lower than
                        the setting will be filtered, default:10
  --min_ava_BQ MIN_AVA_BQ
                        Minimum average base quality required to consider a
                        base for calling, default:20
  --minCov_for_snp MINCOV_FOR_SNP
                        Minimum read counts required to call a snp, default:2
  --minCov_for_indel MINCOV_FOR_INDEL
                        Minimum read counts required to call a indel,
                        default:4
  --minRatio_snp MINRATIO_SNP
                        Minimum variant supported read count ratio for SNP,
                        default:0.125000
  --minRatio_indel MINRATIO_INDEL
                        Minimum variant supported read count ratio for INDEL,
                        default:0.100000
  --min_r_snp_HC MIN_R_SNP_HC
                        Min supporting reads count ratio to detect a high-
                        confidence SNP variant candidate, default:0.300000
  --min_r_indel_HC MIN_R_INDEL_HC
                        Min supporting reads count ratio to detect a high-
                        confidence indel variant candidate, default:0.500000
  --min_c_snp_HC MIN_C_SNP_HC
                        Min supporting reads count to detect a high-confidence
                        SNP variant candidate, default:12
  --min_c_indel_HC MIN_C_INDEL_HC
                        Min supporting reads count to detect a high-confidence
                        indel variant candidate, default:10
  --chrName CHRNAME     The name of reference to be processed, default:22
  --chrStart CHRSTART   The 1-based starting positions of the reference to be
                        processed
  --chrEnd CHREND       The inclusive ending positions of the reference to be
                        processed
  --variantType VARIANTTYPE
                        Extract candidates of SNP, indel or all of small
                        variant,[snp, indel, all], default: all
```

**localMSA.py**

```
  -h, --help            show this help message and exit
  --fin_bam FIN_BAM     Sorted bam file, default: input_local.bam
  --fin_ref FIN_REF     Reference fasta input
  --fout_vcf FOUT_VCF   vcf path for detected variants, default: detect.vcf
  --fin_can FIN_CAN     Candidate generated by ExtractVariantCandidate.py,
                        default: PIPE
  --fin_repeats FIN_REPEATS
                        Tandem repeat intervals for indel assembly, default:
                        None
  --chrName CHRNAME     The name of reference to be processed, default: None
  --chrStart CHRSTART   The 1-based starting positions of the reference to be
                        processed
  --chrEnd CHREND       The inclusive ending positions of the reference to be
                        processed
  --perror_for_snp PERROR_FOR_SNP
                        P-error is the probability of observing a heterozygote
                        SNP, default:0.100000
  --perror_for_indel PERROR_FOR_INDEL
                        P-error is the probability of observing a heterozygote
                        indel, default:0.100000
  --ratio_identity_snp RATIO_IDENTITY_SNP
                        min frequency of each consensus for heterozygote SNP,
                        default:0.200000
  --ratio_identity_indel RATIO_IDENTITY_INDEL
                        min frequency of each consensus for heterozygote
                        indel, default:0.200000
  --minMQ MINMQ         Minimum Mapping Quality. Mapping quality lower than
                        the setting will be filtered, default:10
  --minCNT MINCNT       Minimum read counts required to call a variant,
                        default:3
  --max_merge_dis MAX_MERGE_DIS
                        Max distance to merge two variant candidates,
                        default:10
  --shift SHIFT         The distance bewteen the detect variant and the
                        activate region, default:5
  --flanking FLANKING   Flanking base pairs around variant site, default: 70
  --useBaseQuality      Use base quality to call variant, which will
                        realignmet local reads to haplotypes and cost more
                        time with little improvment of performance, default:
                        True
  --useAllReads         Use all local reads (including spanning reads and
                        part-overlapped reads) to perform MSA, default: True
  -x MISMATCH, --mismatch MISMATCH
                        Mismatch penalty. default: 4
  -d MISMATCH2, --mismatch2 MISMATCH2
                        Mismatch penalty. default: 6
```

## Datasets

We implemented the benchmarks on a server with AMD Ryzen 3950X CPU and 120GB RAM, running Linux Ubuntu 16.04.

The ground truth variants for simulation generated by high confidence variants of HG001 sample from GIAB, which can be downloaded from:
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh37

The real HG002 Illumina PE150 reads was downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/

The real HG002 Illumina PE250 reads was downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps

The GIAB ground truth set and corresponding high confidence region set were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/

## Contact

For advising, bug reporting and requiring help, please post on **[Github Issues](https://github.com/ydLiu-HIT/Psi-caller/issues)** or contact [ydliu@hit.edu.cn](mailto:ydliu@hit.edu.cn).
