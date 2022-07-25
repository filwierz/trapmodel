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
parser.add_argument('--rm', type=argparse.FileType('r'), default=None,dest="rm", required=True, help="RepeatMasker's .out file")
#  SW   perc perc perc  query     position in query              matching     repeat          position in repeat
#score   div. del. ins.  sequence  begin    end          (left)   repeat       class/family begin   end    (left)     ID
#
#71128    2.0  1.3  1.1  2L_RaGOO        10     8592 (24528620) + DM14101      Unspecified     924   9525  (1129)     1
parser.add_argument("--minlen", type=int, required=False, dest="minlen", default=0, help="minimum length of repeatmasker feature")
parser.add_argument("--maxdiv", type=float, required=False, dest="maxdiv", default=100.0, help="maximum divergence of repeatmasker feature")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="minimum mapping quality of small RNA read")
parser.add_argument("--window", type=int, required=False, dest="window", default=0, help="window size of butterfly signal")
parser.add_argument("--id", type=str, required=True, dest="sampleid", default=None, help="id of the sample")


args=parser.parse_args()

ws=args.window
mmq=args.minmq

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
    else:
        toit[chromo][pos][0]+=1


for k in args.rm:
    k=k.rstrip("\n").lstrip(" ")
    r=re.split("\s+",k)
    if r[0]=="SW" or r[0]=="score" or r[0]=="":
        continue
    if float(r[1]) > args.maxdiv:
        continue
    if (int(r[6])-int(r[5])+1) < args.minlen:
        continue

    chro=str(r[4])
    start=int(r[5])
    end=int(r[6])
    TE=str(r[9])

    lw=start-ws
    rw=end+ws


    for c,tmp in toit.items():
        if c==chro:
            for p,ch in tmp.items():
                if lw<=p<=start:
                    topr[TE][chro][start][end][0]+=ch[0]
                    topr[TE][chro][start][end][1]+=ch[1]
                if end<=p<=rw:
                    topr[TE][chro][start][end][2]+=ch[0]
                    topr[TE][chro][start][end][3]+=ch[1]


for t,tmp in topr.items():
    for c,tmp2 in tmp.items():
        for s,tmp3 in tmp2.items():
            for e,ch in tmp3.items():
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(t,c,s,e,ch[0],ch[1],ch[2],ch[3],args.sampleid))
       



    
    