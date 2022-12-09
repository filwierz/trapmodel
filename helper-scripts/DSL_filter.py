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

parser.add_argument('--bed', type=argparse.FileType('r'), default=None,dest="bed", required=False, help="a cluster bed file")
parser.add_argument('--dsl', type=argparse.FileType('r'), default=None,dest="dsl", required=True, help="dsl signatures")
parser.add_argument("--th", type=int, required=False, dest="th", default=0, help="minimum piRNA counts")


args=parser.parse_args()

clu=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str]))))
topr=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[float,float,float,float,str])))))

if args.bed is not None:    
    for l in args.bed:
        #2L_RaGOO        2    19884                 2L_TAS   1000      False
        l=l.rstrip("\n")
        m=l.split("\t")
    
        chro=str(m[0])
        start=int(m[1])+1
        end=int(m[2])+1
        cl=str(m[3])
    
        clu[chro][start][end][0]=cl

for line in args.dsl:
    #DM06920 2L_RaGOO   11685   13214  0.000000000000000 0.000000000000000  0.926303308755419  0.000000000000000 Canton-S
    line=line.rstrip("\n")
    n=line.split("\t")
    te=str(n[0])
    chromo=str(n[1])
    st=int(n[2])
    en=int(n[3])
    ls=float(n[4])
    la=float(n[5])
    rs=float(n[6])
    ra=float(n[7])
    strain=str(n[8])

    if la<args.th or rs<args.th:
        continue

    topr[chromo][st][en][te][0]=ls
    topr[chromo][st][en][te][1]=la
    topr[chromo][st][en][te][2]=rs
    topr[chromo][st][en][te][3]=ra
    topr[chromo][st][en][te][4]=strain

    for c,tmp in clu.items():
        if chromo==c:
            for s,tmp2 in tmp.items():
                if s<=st or s<=en:
                    for e,ch in tmp2.items():
                        if st<=e or en<=e:
                            del topr[chromo][st][en][te]



for c,tmp in topr.items():
    for s,tmp2 in tmp.items():
        for e,tmp3 in tmp2.items():
            for t,ch in tmp3.items():
                print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(t,c,s,e,ch[0],ch[1],ch[2],ch[3],ch[4]))
