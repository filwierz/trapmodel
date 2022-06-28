import os
import sys
import re
import argparse
import random
import inspect
import collections
#import numpy
import math
import multiprocessing

import time

def current_milli_time():
    return round(time.time() * 1000)

# java -jar invade-v801.jar cluster  --x 0.01 --max-ins 10000 --min-w 0.1 --genome mb:10,10,10,10,10
# --rr cm_mb:4,4,4,4,4 --cluster each:kb:100 --u 0.1 --N 1000 --gen 5000 --basepop seg:100 --rep 100
# --silent --steps 20  --no-x-cluins --simid 'kb100'  > diffclusi/kb100.txt &
                                    
def get_basis(invade):
    return "java -Xmx4g -jar {0} cluster --ignore-failed --no-x-cluins --gen 10000 --steps 100 --max-ins 10000 --min-w 0.1 --genome kb:40000,40000,40000,40000,40000 --rr cm_mb:4,4,4,4,4 --cluster each:kb:1400 --ref each:kb:1400 --rep 1 --silent".format(invade)

def get_filter():
    return "| grep -v \"^#\" | awk '$28!=\"base\"' "

def get_rand_u():
    # from 0.005 to 0.5 uniformly log distributed
    #return 10**random.uniform(-0.301029995664,-2.30103)
    return 0.1

#def get_rand_x():
    # from 0.0001 to 0.5 uniformly log distributed
    #return 10**random.uniform(-0.301029995664,-4)
    #return 0.0
    


def run_cluster_negsel(invade,count,output):
    """
    TE invasion that is stopped by cluster insertions and neg selection against TEs
    """
    basis =get_basis(invade) 
    commandlist=[]
    for i in range(0,count):
        #x=get_rand_x()
        u=get_rand_u()
        tr=current_milli_time()+i
        command=basis+" --nsmodel site:0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 --u {0} --replicate-offset {1} --tally-file tally{1}.txt --sfs-file sfs{1}.txt --hohe-file hohe{1}.txt --mhp-file mhp{1}.txt --N 1000 --basepop seg:1000 --seed {2} ".format(u,i,tr)
        ri=random.random()
        command+= "--simid \"{0}\"  {1} > {2}{3}".format(u,get_filter(),output,i)
        commandlist.append(command)
    return commandlist;


parser = argparse.ArgumentParser(description="""           
Description
-----------
    Simulation storm""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Robert Kofler
    Filip Wierzbicki 
""")


parser.add_argument("--number", type=int, required=True, dest="count", default=None, help="the number of simulations")
parser.add_argument("--threads", type=int, required=True, dest="threads", default=None, help="the threads of simulations")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the outputfile of simulations")
parser.add_argument("--invade", type=str, required=True, dest="invade", default=None, help="the invade.jar")
# parser.add_argument("--type", type=str, required=True, dest="type", default=None, help="nostop|ns|nst|cluster|nscluster")
parser.add_argument("--silent",action="store_true", dest="silent", help="be quit; default=False")



args = parser.parse_args()
#0   1    2    3     4     5    6       7    8    9           10     11   1     13   14        15   16   17   18   19   20       21   22        23 
#1	0	|	0.01	1.00	0.01	0.0005	0	|	0.0000	0.00	0.0000	0	|	0.0000	10	0	0.10	0.00	1000	|	u=0.01	ci=1		1
#1	100	|	0.11	1.00	0.12	0.0074	0	|	0.0000	0.00	0.0000	0	|	0.0010	8	0	0.35	0.00	1000	|	u=0.01	ci=1		1
#1	200	|	0.40	1.00	0.50	0.0075	0	|	0.0270	0.03	0.0070	0	|	0.0060	33	2	0.70	0.17	1000	|	u=0.01	ci=1		1
#1	300	|	0.84	1.00	1.76	0.0092	0	|	0.2700	0.30	0.0370	0	|	0.0090	96	4	1.33	0.51	1000	|	u=0.01	ci=1		1

commandlist=run_cluster_negsel(args.invade,args.count,args.output)

"""
if(args.type=="nostop"):
    commandlist=run_nostop(args.invade,args.count,args.output)
elif(args.type=="ns"):
    commandlist=run_negsel(args.invade,args.count,args.output)
elif(args.type=="nst"):
    commandlist=run_negsel_t(args.invade,args.count,args.output)
elif(args.type=="cluster"):
    commandlist=run_cluster(args.invade,args.count,args.output)
elif(args.type=="nscluster"):
    commandlist=run_cluster_negsel(args.invade,args.count,args.output)
else:
    raise Exception("Unknown simulation type")
"""
    
    



 
def submit_job_max_len(commandlist, max_processes):
    import subprocess
    import time
    sleep_time = 10.0
    processes = list()
    for command in commandlist:
        if(not args.silent):
            print ('running {n} processes. Submitting {proc} '.format(n=len(processes),proc=str(command)))
        processes.append(subprocess.Popen(command, shell=True, stdout=None))
        while len(processes) >= max_processes:
            time.sleep(sleep_time)
            processes = [proc for proc in processes if proc.poll() is None]
    while len(processes) > 0:
        time.sleep(sleep_time)
        processes = [proc for proc in processes if proc.poll() is None]
 
submit_job_max_len(commandlist, max_processes=args.threads)
print ("Done")
    


    

