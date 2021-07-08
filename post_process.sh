#!/bin/bash

prefix=$1
output_vcf=$2

for i in ${prefix}/var.*.vcf; do awk '!a[$1" "$2]++{print}' $i > ${prefix}/tmp.vcf ; mv ${prefix}/tmp.vcf ${i}; done 

vcfcat ${prefix}/var.*.vcf | bcftools sort -m 2G >${output_vcf}

#rm ${prefix}/var.*.vcf


