#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 11:04:46 2018

@author: antoine
"""
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import stats as st
import mmap
from tqdm import *



def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def firstselect(liste,thresh,graph):
    cpt=thresh+1
    i=200
    while cpt > thresh:
        x=[]
        cpt=0
        for n,_ in liste:
            if (n,n) in graph.edges:
                if graph.node[n]['outscore']*graph[n][n]['weight']>i and len(graph[n][n]['peaks'])>1 and graph[n][n]['site']!='prom':
                    x.append(n)
                    cpt+=1
        i+=1
    return x,i

def main():

    filep=sys.argv[1]
    #filep="/home/antoine/NAS_Public/data/projects/Antoine_Networks/results/outputDir/ATAC.H1hES.rep_1.bwa.q20.rmdup_peaks.narrowPeak.cross.bed"
    
    print(time.ctime())
    print("building network... 3/8")
    G=nx.DiGraph() 
    TF_binding=set()
    
    
    with open(filep,"r") as cross:
        for line in tqdm(cross, total=get_num_lines(filep)):
#        for line in cross.readlines():
            _,_,_,allgenes,score,_,peak,types = line.strip().split("\t")
            tf=allgenes.split("_")[0]
            TF_binding.add(tf)
            genes=allgenes.split("_")[1].split(",")
            for i in genes:
                if (tf,i) not in G.edges():
                    G.add_edge(tf,i, site=types, peaks=set(), motifs=float(score), weight=float(score))
                    G[tf][i]['peaks'].add(peak)
                else:
                    G[tf][i]['peaks'].add(peak)
                    G[tf][i]['motifs']+=float(score)
                    G[tf][i]['weight']=G[tf][i]['motifs']*len(G[tf][i]['peaks'])
                    if G[tf][i]['site']!=types:
                        G[tf][i]['site']='both'
    
    
    print(time.ctime())
    print("calculating scores... 4/8")
    select=[]
    select2=[]
    links=[]            
    for n in G:                
        G.node[n]['inscore']=float(np.log2(0.001+(sum(G[i][j]['weight'] for i,j in G.in_edges(n)))))
        
#        if ("NANOG",n) in G.edges():
#            tt.append((n,float(np.log2(0.001+(G["NANOG"][n]['weight'])))))
        
        if n in TF_binding:
            #Outscore only on TFs
    #        G.node[n]['outTF']=0
    #        for j in G:
    #            if j in TF_binding and (n,j) in G.edges():
    #                G.node[n]['outTF']+=G[n][j]['weight']
    #        G.node[n]['outTF']=float(np.log2(0.001+G.node[n]['outTF']))
            for j in G:
                if j in TF_binding and (n,j) in G.edges():
                    if (j,n) in G.edges():
                        if (n,j) not in links and (j,n) not in links:
                            links.append((((n,j)),G[n][j]['weight']+G[j][n]['weight']))
                    else:
                        if (n,j) not in links:
                            links.append((((n,j)),G[n][j]['weight']))
    
            G.node[n]['outscore']=float(np.log2(0.001+(sum(G[i][j]['weight'] for i,j in G.out_edges(n)))))
            if (n,n) in G.edges():
                select.append((n,G.node[n]['outscore']))
                if float(np.log2(0.001+(G[n][n]['weight'])))>0:
                    select2.append((n,float(np.log2(0.001+(G[n][n]['weight'])))))
            G.node[n]['IsTF']='yes'
        else:
            G.node[n]['outscore']=0.0
            G.node[n]['IsTF']='no'    
    
             
    selected = sorted(select, key=lambda tup: tup[1])
    tfname, outsc = zip(*selected)
    colors=[]
    place1=[]
    name1=[]
    for i in tfname:
        if i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="KLF7" or i=="ISL1" or i=="PHOX2B": #i=="ZNF281" or i=="NANOG" or i=="POU5F1" or i=="SOX2":
            colors.append('r')
            place1.append(tfname.index(i))
            name1.append(i)
        else:
            colors.append('b')
    indices = np.arange(len(selected))
    plt.figure()
    plt.bar(indices, outsc, color=colors)
    plt.xticks(place1, name1, rotation='vertical')
    plt.tight_layout()
    plt.show(block=False)     
    
    selected2 = sorted(select2, key=lambda tup: tup[1])
    tfname2, autow = zip(*selected2)
    colors2=[]
    place2=[]
    name2=[]
    for i in tfname2:
        if i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="KLF7" or i=="ISL1" or i=="PHOX2B": #i=="ZNF281" or i=="NANOG" or i=="POU5F1" or i=="SOX2":
            colors2.append('r')
            place2.append(tfname2.index(i))
            name2.append(i)
        else:
            colors2.append('b')
    indices2 = np.arange(len(selected2))
    plt.figure()
    plt.bar(indices2, autow, color=colors2)
    plt.xticks(place2, name2, rotation='vertical')
    plt.tight_layout()
    plt.show(block=False)
   
#    links2 = sorted(links, key=lambda tup: tup[1])
#    lname, allinks = zip(*links2)
#    colors=[]
#    place=[]
#    name=[]
#    for i,j in lname:
#        if (i=="NANOG" and (j=="POU5F1" or j=="SOX2")) or (i=="POU5F1" and j=="SOX2"):  #i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="TWIST1" or i=="MYCN" or i=="PHOX2A" or i=="PHOX2B":
#            colors.append('r')
#            place.append(lname.index((i,j)))
#            name.append((i,j))
#        else:
#            colors.append('b')
#    indices = np.arange(len(links2))
#    plt.figure()
#    plt.bar(indices, allinks, color=colors)
#    plt.xticks(place, name, rotation='vertical')
#    plt.tight_layout()
#    plt.show(block=False)
    
    
    print(time.ctime())
    print("finding cliques... 5/8")
    
    seuil=70
    x,i = firstselect(select,seuil,G)
    Gc=nx.Graph()
    for n in x:
        Gc.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        for m in x:
            if (n,m) in G.edges and (m,n) in G.edges:
                Gc.add_edge(n, m, site=G[n][m]['site'], peaks=G[n][m]['peaks'], 
                        motifs=G[n][m]['motifs'], weight=G[n][m]['weight'])
#    Gclique=nx.find_cliques(Gc)
#    a=len(list(Gclique))
#                    
#    while a < 25000:
#        print(a, seuil)
#        if a < 1000:
#            seuil=seuil+10
#        else:
#            seuil=seuil+2
#        x,i = firstselect((x,i),seuil,G)
#        Gc.clear()
#        for n in x:
#            Gc.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
#            for m in x:
#                if (n,m) in G.edges and (m,n) in G.edges:
#                    Gc.add_edge(n, m, site=G[n][m]['site'], peaks=G[n][m]['peaks'], 
#                            motifs=G[n][m]['motifs'], weight=G[n][m]['weight'])
#        print("graph done")
#        Gclique=nx.find_cliques(Gc)
#        a=len(list(Gclique))
    
    print(x)
    print(len(x), i-1)
    
    
    #    #out=[]
    #    names=[]
    #    ins=[]            
    #    #for n in Gc:
    #    #    out.append(Gc.node[n]['outscore'])
    #    for e1, e2 in G.edges:
    #        if G[e1][e2]['weight']>400:
    #            ins.append(G[e1][e2]['weight'])
    #            names.append(e1+"+"+e2)
    #    insf, namesf = zip(*sorted(zip(ins, names)))
    #    #plt.figure()
    #    #plt.hist(out)
    #    plt.figure()
    #    plt.plot(insf)
    #    plt.xticks(range(len(namesf)), namesf, rotation=45)
    #    plt.show()
    
    
    score={}
    clique1=[]
    clique2=[]
    clique3=[]
    maxi1=0
    maxi2=0
    maxi3=0
    i=0
    
    print(time.ctime())
    print("finding core TFs... 6/8")
    
    Gclique=nx.find_cliques(Gc)  
    for n in Gclique:
        i+=1
        tmps=(sum(G[m1][m2]['weight'] for m1 in n for m2 in n))/(len(n)*len(n))
        tmpout=(sum(G.node[m]['outscore'] for m in n))/len(n)
        Cscore=tmps*tmpout
        if Cscore>maxi3:
            if Cscore>maxi2:
                if Cscore>maxi1:
                    maxi3=maxi2
                    clique3=clique2
                    maxi2=maxi1
                    clique2=clique1
                    maxi1=Cscore
                    clique1=n
                else:
                    maxi3=maxi2
                    clique3=clique2
                    maxi2=Cscore
                    clique2=n
            else:
                maxi3=Cscore
                clique3=n
    print("nombre de cliques :", i)        
    
    print(time.ctime())
    print("finding most regulated genes... 7/8")    
    reg={}
    listclique=[clique1,clique2,clique3]
    for cl in listclique:    
        for n in cl:
            p=str(listclique.index(cl)+1)
            reg[p]=reg.get(p, {})
            for m in G.neighbors(n):
                reg[p][m]=reg[p].get(m, {})
                reg[p][m]['liste']=reg[p][m].get('liste', set())
                if (n,m) in G.edges:
                    reg[p][m]['liste'].add(n)
                    if 'score' in reg[p][m]:
                        reg[p][m]['score']+=G[n][m]['weight']
                    else:
                        reg[p][m]['score']=G[n][m]['weight']
    
    Gcl1=nx.DiGraph()
    for n in clique1:
        Gcl1.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique1:
            Gcl1.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            motifs=G[n][m]['motifs'], weight=G[n][m]['weight'])
    
    Gcl2=nx.DiGraph()
    for n in clique2:
        Gcl2.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique2:
            Gcl2.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            motifs=G[n][m]['motifs'], weight=G[n][m]['weight'])
            
    Gcl3=nx.DiGraph()
    for n in clique3:
        Gcl3.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique3:
            Gcl3.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            motifs=G[n][m]['motifs'], weight=G[n][m]['weight'])
    
#    d={}
#    f=open("/home/antoine/NAS_Public/data/projects/Robel_NB/MFischer_tumor_expression_data/MFischer_RNAseq_GSE49711_SEQC_NB_TUC_G_log2.NORM.txt")
#    for line in f.readlines()[1:]:
#        d[line.split("\t")[0]]=d.get(line.split("\t")[0], [])
#        for i in line.split("\t")[1:]:
#            d[line.split("\t")[0]].append(i)
#        d[line.split("\t")[0]]=np.array(d[line.split("\t")[0]])
#    
#    cor={}
#    pairs={}   
#    for n in TF_binding:
#        for m in TF_binding:
#            if m != n and (n,m) in G.edges() and m in d and n in d:
#                if str(m)+"-"+str(n) in pairs:
#                    pairs[str(m)+"-"+str(n)]+=G[n][m]['weight']
#                else:
#                    pairs[str(n)+"-"+str(m)]=G[n][m]['weight']
#                if str(m)+"-"+str(n) not in cor:
#                    cor[str(n)+"-"+str(m)]=abs(st.spearmanr(d[n],d[m])[0])
#    d.clear()
#    na=[]
#    we=[]
#    co=[]
#    for i in pairs:
#        na.append(i)
#        we.append(pairs[i])
#    pw, pn = zip(*sorted(zip(we,na)))
#    for i in pn:
#        co.append(cor[i])
#    co1=co[0:int(len(co)/3)]
#    co2=co[int(len(co)/3):int(2*len(co)/3)]
#    co3=co[int(2*len(co)/3):]
#    print(len(co),len(pn),len(TF_binding))
#    plt.figure()
#    plt.plot(co, ".") 
#    plt.show(block=False)     
    
    
    #    edg=[]
    #    names=[]        
    #    for n in G2:
    #        if (n,n) in G2.edges:
    #            edg.append(G[n][n]['weight'])
    #            names.append(n)
    #    edgf, namesf = zip(*sorted(zip(edg, names)))
    #    plt.figure()
    #    plt.plot(edgf)
    #    plt.xticks(range(len(namesf)), namesf, rotation=45)
    #    plt.show()
      
    print(time.ctime())     
    print("writing files... 8/8")
    
    Gsub_path=filep.replace(".cross.bed",".gene_list")
    for n in reg:
        with open(Gsub_path+"."+str(n)+".txt","w") as file:
           for i in reg[n]:
                file.write(i+"\t"+str(reg[n][i]['score'])+"\t"+G.node[i]['IsTF']+"\t"+','.join(reg[n][i]['liste'])+"\n")
    
    nx.write_gexf(Gcl1,Gsub_path+".1.gexf")
    nx.write_adjlist(Gcl1,Gsub_path+".1.adjlist")
    nx.write_gexf(Gcl2,Gsub_path+".2.gexf")
    nx.write_adjlist(Gcl2,Gsub_path+".2.adjlist")
    nx.write_gexf(Gcl3,Gsub_path+".3.gexf")
    nx.write_adjlist(Gcl3,Gsub_path+".3.adjlist")
    
    print(time.ctime())
    print("done !")
    plt.show()

if __name__ == '__main__':
    main()