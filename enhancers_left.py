#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def main():
    enhancers_f5 = set()
    enhancers_fantom5 = open(sys.argv[1])
    enhancers_all = open(sys.argv[2])
    

    for line in enhancers_fantom5.readlines():
        chrom, start, end, peak_name, gene_name = line.strip().split("\t")
        if peak_name not in enhancers_f5: 
            enhancers_f5.add(peak_name)
    

    for line2 in enhancers_all.readlines():
        chrom, start, end, kind, peak_name_a, mid =  line2.strip().split("\t")
        if peak_name_a not in enhancers_f5:
            print(chrom, start, end, kind, peak_name_a, mid, sep="\t")

if __name__ == '__main__':
    main()
