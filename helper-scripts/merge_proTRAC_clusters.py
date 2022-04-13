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

parser.add_argument('--bed', type=argparse.FileType('r'), default=None,dest="bed", required=True, help="a positionally sorted minimal bed files")
parser.add_argument('--polyn', type=argparse.FileType('r'), default=None,dest="polyn", required=True, help="a polyN file from cusco")

args=parser.parse_args()
newcl=collections.defaultdict((lambda:collections.defaultdict(lambda:[int,str])))

chromosomes=[]

for l in args.bed:
    l=l.rstrip("\n")
    m=l.split("\t")
    
    chro=str(m[0])
    start=int(m[1])
    end=int(m[2])
    cl=str(m[3])
    
    if chro not in chromosomes:
        newcl[chro][start][0]=end
        newcl[chro][start][1]=cl
        
        chromosomes.append(chro)
        
    else:
        for c,tmp in newcl.items():
            if c==chro:
                last=max(tmp, key=tmp.get)
                l1=tmp[last][0]-last
                gap=start-tmp[last][0]
                l2=end-start
                if (l1+l2)>gap:
                    newcl[c][last][0]=end
                    break
                else:
                    newcl[chro][start][0]=end
                    newcl[chro][start][1]=cl
                    break

gapcl=[]
                  
for c,tmp in newcl.items():
    for pos,ch in tmp.items():
        for r in args.polyn:
            r=r.rstrip("\n")
            t=r.split("\t")
            chrom=str(t[0])
            gs=int(t[1])
            ge=int(t[2])
                   
            if chrom==c and (pos<=(gs or ge)<=ch[0]):
                gapcl.append(ch[1])         
        args.polyn.seek(0)
        
for c,tmp in newcl.items():
    for pos,ch in tmp.items():
        if ch[1] in gapcl:
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(c,pos,ch[0],ch[1],"0","False"))
        else:
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(c,pos,ch[0],ch[1],"1000","False"))
    