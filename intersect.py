#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

def main():
    open_region = {}
    open_peak = {}
    prev_chrom = None

    for line in sys.stdin.readlines():
        pieces = line.strip().split("\t")
        (chrom, start, end, kind, name), rest = pieces[:5], pieces[5:]
        start = int(start)
        end = int(end)

        if chrom != prev_chrom:
            open_region = {}
            open_peak = {}
        else:
            open_region = {region: region_end for region, region_end in open_region.items() if start <= region_end }

            open_peak = {peak: (peak_start, peak_end) for peak, (peak_start, peak_end) in open_peak.items() if start <= peak_end }


        if kind == sys.argv[1]:
            region = name
            open_region[region] = end
            for peak, (peak_start, peak_end) in open_peak.items():
                print(chrom, peak_start, peak_end, peak, region, sep="\t")

        if kind == 'peak':
            peak, peak_start, peak_end = name, start, end
            open_peak[peak] = peak_start, peak_end
            for region in open_region:
                print(chrom, peak_start, peak_end, peak, region, sep="\t")
        

        prev_chrom = chrom

if __name__ == '__main__':
    main()

