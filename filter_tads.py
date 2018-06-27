#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def main():
    expressed_genes = open(sys.argv[1])
    enhancers_tad = open(sys.argv[2])
    expr_genes = set()
    gene_info = {}
    
    for line in expressed_genes.readlines():
        _,_,_,_, gene_name = line.strip().split("\t")
        for chunk in gene_name.split(";"):
            key, value = chunk.split("=")
            gene_info[key] = value
        expr_genes.add(gene_info["gene_name"])
    
    for line in enhancers_tad.readlines():
        chrom, start, end, peak, genes = line.strip().split("\t")
        gene_names = set(genes.split("|"))
        matching_genes = expr_genes.intersection(gene_names)
        if len(matching_genes) > 0:
            print(chrom, start, end, peak, "|".join(matching_genes), sep="\t")

if __name__ == '__main__':
    main()

