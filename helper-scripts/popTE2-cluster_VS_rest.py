from __future__ import division

#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import collections


from collections import defaultdict

parser = argparse.ArgumentParser(description="""           
Description
-----------
    test""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 2.7

""")


parser.add_argument('--popte2', type=argparse.FileType('r'), default=None,dest="pt", required=True, help="popte2 te insertion output")
parser.add_argument('--pic', type=argparse.FileType('r'), default=None,dest="pic", required=True, help="4 column 1-based piRNA cluster annotations of a genome")
parser.add_argument('--minfreq', type=float, required=False, dest="minfreq", default=0.0, help="min freq of TE insertion")

args = parser.parse_args()
fre=args.minfreq
ptopr = collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str,str,str,str,float,str,str])))))


for line in args.pt:
    #1       arm_2L  7101362 .       I-element       non-LTR FR      -       0.951
    line=line.rstrip("\n")
    a=line.split("\t")
    
    
    for row in args.pic:
        #tig00004115_pilon_pilon 1553768 1832981 cl1 
        row=row.rstrip("\n")
        b=row.split(" ")

        if (str(a[1])==str(b[0]) and int(b[1])<=int(a[2]) and int(a[2])<=int(b[2])):
            if (float(a[8])>=fre):
                
                ptopr[a[4]][a[1]][a[2]][a[0]][0]=a[3]
                ptopr[a[4]][a[1]][a[2]][a[0]][1]=a[5]
                ptopr[a[4]][a[1]][a[2]][a[0]][2]=a[6]
                ptopr[a[4]][a[1]][a[2]][a[0]][3]=a[7]
                ptopr[a[4]][a[1]][a[2]][a[0]][4]=a[8]
                ptopr[a[4]][a[1]][a[2]][a[0]][5]="cl"
                ptopr[a[4]][a[1]][a[2]][a[0]][6]=b[3]
                break
        else:
            if (float(a[8])>=fre):
                
                ptopr[a[4]][a[1]][a[2]][a[0]][0]=a[3]
                ptopr[a[4]][a[1]][a[2]][a[0]][1]=a[5]
                ptopr[a[4]][a[1]][a[2]][a[0]][2]=a[6]
                ptopr[a[4]][a[1]][a[2]][a[0]][3]=a[7]
                ptopr[a[4]][a[1]][a[2]][a[0]][4]=a[8]
                ptopr[a[4]][a[1]][a[2]][a[0]][5]="non-cl"
                ptopr[a[4]][a[1]][a[2]][a[0]][6]="-"
                
            
    args.pic.seek(0)


for t,tmp in ptopr.items():
    for cont,tmp2 in tmp.items():
        for st,tmp3 in tmp2.items():
            for sample,ch in tmp3.items():
                print("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(sample,cont,st,ch[0],t,ch[1],ch[2],ch[3],ch[4],ch[5],ch[6]))
        

