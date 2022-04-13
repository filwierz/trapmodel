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

parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="a sam file from bwasw")
parser.add_argument('--polyn', type=argparse.FileType('r'), default=None,dest="polyn", required=True, help="a polyN file from cusco")
parser.add_argument('--fai', type=argparse.FileType('r'), default=None,dest="fai", required=True, help="a fasta index file from samtools")
parser.add_argument('--mq', type=int, default=2,dest="mqt", required=False, help="minimum mapping quality of flanking sequence")

args=parser.parse_args()
newcl=collections.defaultdict((lambda:[int,int,int,int,str]))

for l in args.sam:
    l=l.rstrip("\n")
    m=l.split("\t")
    if m[0]=="@SQ":
        continue
    if m[2]=="*":
        continue
    else:
        flank=str(m[0])
        flag=int(m[1])
        chro=str(m[2])
        pos=int(m[3])
        mq=int(m[4])
        leng=len(str(m[9]))
        if mq<args.mqt:
            continue
        if flank not in newcl.keys():
            newcl[flank][0]=flag
            newcl[flank][1]=chro
            newcl[flank][2]=pos
            newcl[flank][3]=mq
            newcl[flank][4]=leng
        elif mq>newcl[flank][3]:
            newcl[flank][0]=flag
            newcl[flank][1]=chro
            newcl[flank][2]=pos
            newcl[flank][3]=mq
            newcl[flank][4]=leng
        else:
            continue
         
for c,ch in newcl.items():
    if c=="X_TAS" or c=="2L_TAS" or c=="3L_TAS":
        if ch[0]==0:
            start=1
            for r in args.polyn:
                    r=r.rstrip("\n")
                    t=r.split("\t")
                    chro=str(t[0])
                    end=int(t[2])
                    
                    #if chro==ch[1] and end<ch[2] and start>1:
                       # start=end
                    if chro==ch[1] and end<ch[2] and end>start:
                        start=end
            args.polyn.seek(0)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(ch[1],start,ch[2],c,"1000","False"))
            
        elif ch[0]==16:
            for i in args.fai:
                i=i.rstrip("\n")
                f=i.split("\t")
                chromo=str(f[0])
                if chromo==ch[1]:
                    end=int(f[1])
                    
            args.fai.seek(0)
            for r in args.polyn:
                    r=r.rstrip("\n")
                    t=r.split("\t")
                    chro=str(t[0])
                    sta=int(t[1])

                    #if chro==ch[1] and sta>ch[2] and end is None:
                       # end=sta
                    if chro==ch[1] and sta>ch[2] and sta<end:
                        end=sta
            args.polyn.seek(0)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(ch[1],ch[2]+ch[4],end,c,"1000","True"))
        
        else:
            raise Exception("unexpected sam flag")
        
    
        
    if c=="2R_TAS" or c=="3R_TAS":
        
        if ch[0]==0:
            start=ch[2]+ch[4]
            
            for i in args.fai:
                i=i.rstrip("\n")
                f=i.split("\t")
                chromo=str(f[0])
                if chromo==ch[1]:
                    end=int(f[1])
                    
            args.fai.seek(0)
            
            for r in args.polyn:
                    r=r.rstrip("\n")
                    t=r.split("\t")
                    chro=str(t[0])
                    sta=int(t[1])
                    
                    #if chro==ch[1] and sta>start and end is None:
                       # end=sta
                    if chro==ch[1] and sta>start and sta<end:
                        end=sta
            args.polyn.seek(0)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(ch[1],start,end,c,"1000","False"))
            
        elif ch[0]==16:
            start=1
            for r in args.polyn:
                    r=r.rstrip("\n")
                    t=r.split("\t")
                    chro=str(t[0])
                    end=int(t[2])
                    if chro==ch[1] and end<ch[2] and end>start:
                        start=end
            args.polyn.seek(0)
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(ch[1],start,ch[2],c,"1000","True"))
        
        else:
            raise Exception("unexpected sam flag")
    
                     
                                        
                            