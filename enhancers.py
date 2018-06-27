#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def main():
    promoter_peaks = set()
    promoters = open(sys.argv[1])
    all_peaks = open(sys.argv[2])


    for line in promoters.readlines():
        chrom, start, end, peak_name, gene_name = line.strip().split("\t")
        if peak_name not in promoter_peaks:
            promoter_peaks.add(peak_name)

    for line in all_peaks.readlines():
        pieces = line.strip().split("\t")
        chrom, start, end, kind, peak, mid = pieces
        if peak not in promoter_peaks:
            print(chrom, start, end, kind, peak, mid, sep="\t")
        
if __name__ == '__main__':
    main()	
	
