#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 14:37:23 2018

@author: antoine
"""
from tqdm import *
import sys
import time
import mmap
import matplotlib.pyplot as plt


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

#filep=sys.argv[1]
filep="/home/antoine/NAS_Public/data/projects/Antoine_Networks/results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.summary.bed"

interact=set()
PPI={}
TF_binding=set()

with open(filep,"r") as summ:
    for line in tqdm(summ, total=get_num_lines(filep)):
#        for line in cross.readlines():
        _,_,_,allgenes,_,_,_,_ = line.strip().split("\t")
        tf=allgenes.split("_")[0]
        TF_binding.add(tf)

with open("/home/antoine/NAS_Public/data/projects/Antoine_Networks/data/PPI.txt") as names:
    for line in names.readlines()[1:]:
        _,ensembl,gene_name = line.strip().split(" ")
        ens=ensembl.split("\"")[1]
        gen=gene_name.split("\"")[1]
        if gen in TF_binding:
            PPI[ens]=gen
        
explist=[]
with open("/home/antoine/NAS_Public/data/projects/Antoine_Networks/data/9606.protein.links.full.v10.5.txt") as PP:
    for line in PP.readlines()[1:]:
        pep1,pep2,_,_,_,_,_,_,_,exp,_,_,_,_,_,_ = line.strip().split(" ")
        if int(exp)>=150:
            explist.append(int(exp))
            p1=pep1.split(".")[1]
            p2=pep2.split(".")[1]
            if p1 in PPI and p2 in PPI and p1!=p2:
                interact.add((PPI[p1],PPI[p2]))
                
plt.figure()
plt.hist(explist)
plt.show(block=False)
    
with open("/home/antoine/NAS_Public/data/projects/Antoine_Networks/data/9606.protein.links.full.v10.5.filtered.txt","w") as file:        
    for (i,j) in interact:
        file.write(i+" "+j+"\n")
file.close()