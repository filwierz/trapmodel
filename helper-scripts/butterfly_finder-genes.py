#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import collections

parser = argparse.ArgumentParser(description="""
Description
-----------
    test""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

""")

parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="a cluster bed file")
#HISEQ:298:CALHPANXX:8:1107:4304:84076/1	0	2L	8	0	23M	*	0	0	GCACGACAGAGGAAGCAGAACAG	BBBBBFFFFFFFFFFFFFFFFFF	MD:Z:23	PG:Z:novoalign.1	IH:i:1	NH:i:55	HI:i:1	NM:i:0	UQ:i:0	AS:i:0	ZS:Z:R
#HISEQ:298:CALHPANXX:8:2114:13501:60563/1	0	2L	8	0	23M	*	0	0	GCACGACAGAGGAAGCAGAACAG	BBBBBFFFFFFFFFFFFFFFFFF	MD:Z:23	PG:Z:novoalign	IH:i:1	NH:i:55	HI:i:1	NM:i:0	UQ:i:0	AS:i:0	ZS:Z:R
parser.add_argument('--bed', type=argparse.FileType('r'), default=None,dest="bed", required=True, help="a bed file for gene annotation")
# 4	720796	776070	20at7147
#2R	8864952	8888496	29at7147
parser.add_argument("--minlen", type=int, required=False, dest="minlen", default=0, help="minimum length of repeatmasker feature")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="minimum mapping quality of small RNA read")
parser.add_argument("--window", type=int, required=False, dest="window", default=0, help="window size of butterfly signal")
parser.add_argument("--id", type=str, required=True, dest="sampleid", default=None, help="id of the sample")


args=parser.parse_args()

ws=args.window
mmq=args.minmq
pp=0

topr=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[0,0,0,0])))))
toit=collections.defaultdict((lambda:collections.defaultdict(lambda:[0,0])))


for l in args.sam:
    l=l.rstrip("\n")
    a=l.split("\t")
    if "@" in a[0]:
        continue
    chromo=str(a[2])
    pos=int(a[3])
    mq=int(a[4])
    rl=len(str(a[9]))
    flag=int(a[1])
    if mq<mmq:
        continue
    if rl<23 or rl>29:
        continue
    if flag & 0x10: #rc
        toit[chromo][pos][1]+=1
        pp+=1
    else:
        toit[chromo][pos][0]+=1
        pp+=1


for k in args.bed:
    k=k.rstrip("\n")
    r=k.split("\t")
    
    chro=str(r[0])
    start=int(r[1])
    end=int(r[2])
    gene=str(r[3])
    if (end-start) < args.minlen:
        continue
    start=start+1
    end=end+1

    lw=start-ws
    rw=end+ws


    for c,tmp in toit.items():
        if c==chro:
            for p,ch in tmp.items():
                if lw<=p<=start:
                    topr[gene][chro][start][end][0]+=ch[0]
                    topr[gene][chro][start][end][1]+=ch[1]
                if end<=p<=rw:
                    topr[gene][chro][start][end][2]+=ch[0]
                    topr[gene][chro][start][end][3]+=ch[1]


for t,tmp in topr.items():
    for c,tmp2 in tmp.items():
        for s,tmp3 in tmp2.items():
            for e,ch in tmp3.items():
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(t,c,s,e,1000000*ch[0]/pp,1000000*ch[1]/pp,1000000*ch[2]/pp,1000000*ch[3]/pp,args.sampleid))
       



    
    