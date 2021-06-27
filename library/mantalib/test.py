#!/usr/bin/env python
# coding=utf-8
import sys
from mantaAss_module import mantaAss
from ksw2_module import ksw2_aligner as aligner

fin = sys.argv[1]
fint = sys.argv[2]
minW = int(sys.argv[3])
maxW = int(sys.argv[4])

Reads = list()
with open(fin, 'r') as f:
    for line in f:
        if line.startswith('>'): continue 

        s = bytes(line.strip(), encoding = "utf-8")
        Reads.append(s)
rdl = [len(s) for s in Reads]
print(max(rdl), min(rdl))

target = str()
with open(fint, 'r') as f:
    for line in f:
        if line.startswith('>'): continue

        target = line.strip()

#res = mantaAss(Reads, min(rdl) - 5, max(rdl) + 5)

res = mantaAss(Reads, minW, maxW)
print(res)

for i in res[1]:
    ali = aligner(bytes.decode(i), target)

    print(ali)


