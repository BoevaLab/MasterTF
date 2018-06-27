#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def main():
    TADS_enhancers = open(sys.argv[1])
    prev_key = None
    peak_genes = set()

    for line in TADS_enhancers.readlines():
        chrom, start, end, peak, genes = line.strip().split("\t")
        genes = genes.split("|")
        
        
        key = chrom, start, end, peak
 
        if prev_key and key != prev_key:
            chrom, start, end, peak = prev_key
            print(chrom, start, end, peak, "|".join(peak_genes), sep="\t")
            peak_genes = set() 
        
        for gene in set(genes):
            peak_genes.add(gene)
        
        prev_key = key
         
           
    if prev_key:        
        chrom, start, end, peak = prev_key
        print(chrom, start, end, peak, "|".join(peak_genes), sep="\t")



if __name__ == '__main__':
    main()
