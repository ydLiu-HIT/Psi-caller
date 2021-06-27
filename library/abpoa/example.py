#!/usr/bin/env python
# coding=utf-8
from poa_module import POA
import sys

is_msa = int(sys.argv[1])

a = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
c = b"ACGTACGTACGTATGTACGTACGTACGTACGTACGTACGTACGTACGT"
d = b"ACGTACGTACGTAGGTACGTACGTACGTACGTAAAACGTACGTACGTACGT"

A = [a,a,a,a,a,a,c,c,c,c,c,d]

consensus = POA(A, 1, 0.3, is_msa)

print(consensus)

#for i in consensus:
#    print(i)
