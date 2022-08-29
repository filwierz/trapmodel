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
    """,formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

""")

parser.add_argument('--bed', type=argparse.FileType('r'), default=None,dest="bed", required=True, help="a cluster bed file")
#2L_RaGOO	10	63820	1	1000	False

args=parser.parse_args()
topr=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str,str,str]))))


for l in args.bed:
    l=l.rstrip("\n")
    m=l.split("\t")
    chr=m[0]
    start=m[1]
    end=m[2]
    cl=m[3]
    gap=m[4]
    irc=m[5]
    if "," not in gap:
        topr[chr][start][end][0]=cl
        topr[chr][start][end][1]=gap
        topr[chr][start][end][2]=irc
    if "," in gap:
        gapm=gap.replace(",", " ") 
        gapm=gapm.replace("1000", "ungapped")
        if "0" in gapm:
            topr[chr][start][end][0]=cl
            topr[chr][start][end][1]="0"
            topr[chr][start][end][2]=irc
        else:
            topr[chr][start][end][0]=cl
            topr[chr][start][end][1]="1000"
            topr[chr][start][end][2]=irc

for c,tmp in topr.items():
    for st,tmp2 in tmp.items():
        for en,ch in tmp2.items():
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(c,st,en,ch[0],ch[1],ch[2]))
    

                
