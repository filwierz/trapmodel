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

parser.add_argument('--clu', type=argparse.FileType('r'), default=None,dest="bed", required=True, help="a cluster bed file")
#2L_RaGOO	10	63820	1	1000	False
parser.add_argument('--rm', type=argparse.FileType('r'), default=None,dest="rm", required=True, help="RepeatMasker's .out file")
#  SW   perc perc perc  query     position in query              matching     repeat          position in repeat
#score   div. del. ins.  sequence  begin    end          (left)   repeat       class/family begin   end    (left)     ID
#
#71128    2.0  1.3  1.1  2L_RaGOO        10     8592 (24528620) + DM14101      Unspecified     924   9525  (1129)     1

parser.add_argument('--ref', type=argparse.FileType('r'), default=None,dest="refbed", required=False, help="a reference bed file")
#2R_RaGOO	10	63820	r1	1000	False

parser.add_argument("--minlen", type=int, required=False, dest="minlen", default=0, help="minimum length of repeatmasker feature")
parser.add_argument("--maxdiv", type=float, required=False, dest="maxdiv", default=100.0, help="maximum divergence of repeatmasker feature")

parser.add_argument("--output", type=str, required=False, dest="output", default="", help="output directory")
parser.add_argument("--sample", type=str, required=True, dest="sample", default="", help="sample name")
parser.add_argument("--approach", type=str, required=True, dest="approach", default="", help="name of approach")




args=parser.parse_args()

outputdir=args.output
sample=args.sample
approach=args.approach

toit=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str,int]))))
rest=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str,int]))))
topr=collections.defaultdict((lambda:collections.defaultdict(lambda:[str])))
toprrest=collections.defaultdict((lambda:collections.defaultdict(lambda:[str])))

abupr=collections.defaultdict((lambda:[0]))
aburest=collections.defaultdict((lambda:[0]))



toitref=collections.defaultdict((lambda:collections.defaultdict(lambda:collections.defaultdict(lambda:[str,int]))))
reftopr=collections.defaultdict((lambda:collections.defaultdict(lambda:[str])))
refabupr=collections.defaultdict((lambda:[0]))


for k in args.rm:
    k=k.rstrip("\n").lstrip(" ")
    r=re.split("\s+",k)
    if r[0]=="SW" or r[0]=="score" or r[0]=="":
        continue
    if float(r[1]) > args.maxdiv:
        continue
    if (int(r[6])-int(r[5])+1) < args.minlen:
        continue
    if r[8]=="+":
        remr=r[13].replace('(', '')
        remr=int(remr.replace(')', ''))
        if ((int(r[11])-1) >= args.minlen) or (remr >= args.minlen):
            continue
    if r[8]=="C":
        remr=r[11].replace('(', '')
        remr=int(remr.replace(')', ''))
        if((int(r[13])-1) >= args.minlen) or (remr >= args.minlen):
            continue

    toit[r[4]][r[5]][r[6]][0]=k
    toit[r[4]][r[5]][r[6]][1]=r[14]
    
    rest[r[4]][r[5]][r[6]][0]=k
    rest[r[4]][r[5]][r[6]][1]=r[14]
    
    toitref[r[4]][r[5]][r[6]][0]=k
    toitref[r[4]][r[5]][r[6]][1]=r[14]
    
for l in args.bed:
    l=l.rstrip("\n")
    m=l.split("\t")
    
    chro=str(m[0])
    start=int(m[1])+1
    end=int(m[2])+1
    cl=str(m[3])
    
    for c,tmp in toit.items():
        if chro==c:
            for st,tmp2 in tmp.items():
                for en,ch in tmp2.items():
                    if (start <= int(st) <= end) or (start <= int(en) <= end):
                        topr[cl][ch[1]][0]=ch[0]#first key: cl -for separate cluster outputs; "cluster" -for single
                        del rest[c][st][en]
                        del toitref[c][st][en]
                        
                        r=re.split("\s+",ch[0])
                        abupr[r[9]][0]+=1
                        
                        
if args.refbed is not None:
    
    for l in args.refbed:
        l=l.rstrip("\n")
        m=l.split(" ")
    
        chro=str(m[0])
        start=int(m[1])+1
        end=int(m[2])+1
        cl=str(m[3])
    
        for c,tmp in toitref.items():
            if chro==c:
                for st,tmp2 in tmp.items():
                    for en,ch in tmp2.items():
                        if (start <= int(st) <= end) or (start <= int(en) <= end):
                            reftopr[cl][ch[1]][0]=ch[0]
                            del rest[c][st][en]
                        
                            r=re.split("\s+",ch[0])
                            refabupr[r[9]][0]+=1


                        
for c,tmp in rest.items():
    for st,tmp2 in tmp.items():
        for en,ch in tmp2.items():
            toprrest["rest"][ch[1]][0]=ch[0]
            
            r=re.split("\s+",ch[0])
            aburest[r[9]][0]+=1
                        
                        
for c,tmp in topr.items():
    outname=outputdir+c+".fasta.out"
    output_file = open(outname, 'w')
    for i,ch in tmp.items():
        output_file.write(re.sub('\s+','\t',ch[0])+"\n")
    output_file.close()


for c,tmp in toprrest.items():
    outname=outputdir+c+".fasta.out"
    output_file = open(outname, 'w')
    for i,ch in tmp.items():
        output_file.write(re.sub('\s+','\t',ch[0])+"\n")
    output_file.close()
    
    
    
    
for c,tmp in reftopr.items():
    outname=outputdir+"ref"+c+".fasta.out"
    output_file = open(outname, 'w')
    for i,ch in tmp.items():
        output_file.write(re.sub('\s+','\t',ch[0])+"\n")
    output_file.close()

    
    
outname=sample+"_"+approach+"_summary.forR"
output_file = open(outname, 'w')
for t,ch in abupr.items():
    output_file.write(str(ch[0])+"\t"+sample+"\t"+t+"\t"+"cluster"+"\n")
for t,ch in aburest.items():
    output_file.write(str(ch[0])+"\t"+sample+"\t"+t+"\t"+"non-cluster"+"\n")
for t,ch in refabupr.items():
    output_file.write(str(ch[0])+"\t"+sample+"\t"+t+"\t"+"ref"+"\n")
output_file.close()




    
