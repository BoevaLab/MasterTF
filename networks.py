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
import math
from scipy import stats as st
import mmap
from tqdm import *


### FOR TQDM (Progress bar when reading large files)
def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

### TO SELECT BEST TFs BEFORE CLIQUE ANALYSIS
def firstselect(thresh,graph):
    cpt=thresh+1
    i=0
    while cpt > thresh:
        x=[]
        cpt=0
        for n in graph:
            if (n,n) in graph.edges:
                if graph.node[n]['outscore']*graph[n][n]['weight']>i and len(graph[n][n]['peaks'])>1 and graph[n][n]['site']!='prom':
                    x.append(n)
                    cpt+=1
        i+=0.1
    return x

def main():

    filep=sys.argv[1]       #summary file
    interaction=sys.argv[2]  #True or False
    imgfolder=sys.argv[3]
#    filep="/home/antoine/NAS_Public/data/projects/Antoine_Networks/results/outputDir/DNAse.H1hES.ENCODE_peaks.narrowPeak.summary.bed"
#    interaction="False"
#    ppifile="/".join(filep.split("/")[0:-1]).replace("results/outputDir","data")+("/9606.protein.links.full.v10.5.filtered.txt")
    filec=filep.replace("summary.bed","fprepared.bed")

    corepeaks=set()
    namefile=filep.split("/")[-1].replace("_peaks.narrowPeak.summary.bed","")
    nameppi='noPPI'   #name used for the png files
    edgethresh=0
    
    maxcliques=35000  #Thresholds that heavily determine the pipeline running time
    maxselec=140
    
    if interaction=='True':  #argv[2], whether or not PPI are used
        edgethresh=1
        ppifile=sys.argv[4]
        nameppi='PPI150'
        interact={}
            
        with open(ppifile) as PP: #Protein-Protein Interactions dictionary
            for line in tqdm(PP, total=get_num_lines(ppifile)):
                pep1,pep2 = line.strip().split(" ")
                interact[pep1]=interact.get(pep1, set())
                interact[pep1].add(pep2)
    
    print(time.ctime())
    
    with open(filec,"r") as core:   #list of peaks found in COREs 
        for line in tqdm(core, total=get_num_lines(filec)):
            _,_,_,_,peak,_ = line.strip().split("\t")
            corepeaks.add(peak)
    
    ###BUILDING A NETWORKX OBJECT FROM THE SUMMARY FILE
    
    print("building network... 3/8")
    G=nx.DiGraph() 
    TF_binding=set()
    oldpeak=0
    oldgenes=[]
    curchrom="chrZ"
    
    with open(filep,"r") as summ:
        for line in tqdm(summ, total=get_num_lines(filep)):
            chrom,start,end,allgenes,score,_,peak,types = line.strip().split("\t")
            tf=allgenes.split("_")[0]
            TF_binding.add(tf)
            genes=allgenes.split("_")[1].split(",")
            for i in genes:
                if tf==i:       #MOTIF ANALYSIS IN ALL PEAKS FOR AUTOLINKS
                    if (tf,i) not in G.edges():
                        G.add_edge(tf,i, site=types, peaks=set(), weight=float(score))
                        G[tf][i]['peaks'].add(peak)
                    else:
                        G[tf][i]['peaks'].add(peak)
                        G[tf][i]['weight']+=float(score)
                        if G[tf][i]['site']!=types:
                            G[tf][i]['site']='both'
                    
                    ### PROTEIN-PROTEIN INTERACTION
                    if interaction=="True":  #Motifs scores of TFs are added to all TFs interaction partners
                        if chrom!=curchrom:
                            overlap=[]
                            listTF=[]
                            listscores=[]
                            overlap.append((int(start),int(end)))
                            listTF.append(tf)
                            listscores.append(float(score))
                            focus=0
                            curchrom=chrom
                        else:
                        #add motifs to 'overlap', compare one motif ('focus') to all others in 'overlap'
                        #for each motif, add its score to every TF that can interact with the motif's TF
                        #so iterate on 'focus', add its score to every TF that interact but is not in 'overlap'
                        #and remove from 'overlap' all motifs that don't overlap with 'focus' 
                            overlap.append((int(start),int(end)))
                            listTF.append(tf)
                            listscores.append(float(score))
                            p=0
                            while(p<focus):
                                if(overlap[focus][0]>overlap[p][1]):
                                    overlap.pop(p)
                                    listTF.pop(p)
                                    listscores.pop(p)
                                    focus=focus-1
                                else:
                                    p=p+1
                            while int(start)>overlap[focus][1]:
                                if listTF[focus] in interact:
                                    for p in interact[listTF[focus]]:
                                        if p not in listTF[:-1]:
                                            for i in oldgenes:
                                                if (p,i) not in G.edges():
                                                    G.add_edge(p,i, site=types, peaks=set(), weight=listscores[focus])
                                                    G[p][i]['peaks'].add(oldpeak)
                                                else:
                                                    G[p][i]['peaks'].add(oldpeak)
                                                    G[p][i]['weight']+=listscores[focus]
                                                    if G[p][i]['site']!=types:
                                                        G[p][i]['site']='both'
                                focus=focus+1
                    oldpeak=peak
                    oldgenes=genes
                            
                            
                else:       #MOTIF ANALYSIS ONLY IN PEAKS IN CREAM-CALLED CORES FOR ALL THE OTHER LINKS
                    if peak in corepeaks:
                        if (tf,i) not in G.edges():
                            G.add_edge(tf,i, site=types, peaks=set(), weight=float(score)) 
                            G[tf][i]['peaks'].add(peak)
                        else:
                            G[tf][i]['peaks'].add(peak)
                            G[tf][i]['weight']+=float(score)
                            if G[tf][i]['site']!=types:
                                G[tf][i]['site']='both'
                        
            ### PROTEIN-PROTEIN INTERACTION
                        if interaction=="True": #Motifs scores of TFs are added to all TFs interaction partners
                            if chrom!=curchrom:
                                overlap=[]
                                listTF=[]
                                listscores=[]
                                overlap.append((int(start),int(end)))
                                listTF.append(tf)
                                listscores.append(float(score))
                                focus=0
                                curchrom=chrom
                            else:
                            #add motifs to 'overlap', compare one motif ('focus') to all others in 'overlap'
                            #for each motif, add its score to every TF that can interact with the motif's TF
                            #so iterate on 'focus', add its score to every TF that interact but is not in 'overlap'
                            #and remove from 'overlap' all motifs that don't overlap with 'focus' 
                                overlap.append((int(start),int(end)))
                                listTF.append(tf)
                                listscores.append(float(score))
                                p=0
                                while(p<focus):
                                    if(overlap[focus][0]>overlap[p][1]):
                                        overlap.pop(p)
                                        listTF.pop(p)
                                        listscores.pop(p)
                                        focus=focus-1
                                    else:
                                        p=p+1
                                while int(start)>overlap[focus][1]:
                                    if listTF[focus] in interact:
                                        for p in interact[listTF[focus]]:
                                            if p not in listTF[:-1]:
                                                for i in oldgenes:
                                                    if (p,i) not in G.edges():
                                                        G.add_edge(p,i, site=types, peaks=set(), weight=listscores[focus])
                                                        G[p][i]['peaks'].add(oldpeak)
                                                    else:
                                                        G[p][i]['peaks'].add(oldpeak)
                                                        G[p][i]['weight']+=listscores[focus]
                                                        if G[p][i]['site']!=types:
                                                            G[p][i]['site']='both'
                                    focus=focus+1
                        oldpeak=peak
                        oldgenes=genes
                 
    ###REMOVING LINKS WITH ONLY FEW BINDING SITES AND LOW SCORE
    for i in G:
        if i in TF_binding:
            for j in G:
                if (i,j) in G.edges():
                    if G[i][j]['weight']<=edgethresh:
                        G.remove_edge(i,j)
       
    
    print(time.ctime())
    print("calculating scores... 4/8")
    ### Nodes Scores calculated from Edges Scores
    select=[]
    select2=[]
    links=[]            
    for n in G:                
        G.node[n]['inscore']=float(np.log2(0.001+(sum(G[i][j]['weight'] for i,j in G.in_edges(n)))))
        if n in TF_binding:
            G.node[n]['outscore']=float(np.log2(0.001+(sum(G[i][j]['weight'] for i,j in G.out_edges(n)))))
            G.node[n]['IsTF']='yes'
        else:
            G.node[n]['outscore']=0.0
            G.node[n]['IsTF']='no'    
    
    
    print(time.ctime())
    print("finding cliques... 5/8")
    
    ### First selection so we don't find too many cliques (time consuming)
    ### Probably could be improved
    
    seuil=70        #really stringent start
    x = firstselect(seuil,G)
    
    ### Undirected graph after first selection required for clique discovery
    
    Gc=nx.Graph()
    for n in x:
        Gc.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        for m in x:
            if (n,m) in G.edges and (m,n) in G.edges:
                Gc.add_edge(n, m, site=G[n][m]['site'], peaks=G[n][m]['peaks'], 
                        weight=G[n][m]['weight'])
    
    Gclique=nx.find_cliques(Gc)  #clique discovery
    a=len(list(Gclique))            
        
    #While we don't get 'too many' cliques, do the first selection again, but be less stringent
    while a < maxcliques and seuil < maxselec:  #Thresholds that heavily determine the pipeline running time
        print(a, seuil)
        if a < maxcliques/10:
            seuil=seuil+10
        else:
            seuil=seuil+2
        x = firstselect(seuil,G)
        Gc.clear()
        for n in x:
            Gc.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
            for m in x:
                if (n,m) in G.edges and (m,n) in G.edges:
                    Gc.add_edge(n, m, site=G[n][m]['site'], peaks=G[n][m]['peaks'], 
                            weight=G[n][m]['weight'])
        print("graph done")
        Gclique=nx.find_cliques(Gc)
        a=len(list(Gclique))
    
    print(a, seuil)
    
    
    ### Directed graph for plots
    Ggraph=nx.DiGraph()
    for n in x:
       Ggraph.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
       for m in x:
           if (n,m) in G.edges():
               Ggraph.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            weight=G[n][m]['weight'])
            
    
    ### SCORES DISTRIBUTIONS
    for n in Ggraph:                
        if n in TF_binding:
            for j in Ggraph:
                if (n,j) in Ggraph.edges():
                    links.append((((n,j)),float(np.log2(0.001+(Ggraph[n][j]['weight']/len(Ggraph))))))
    
            if (n,n) in Ggraph.edges():
                select.append((n,Ggraph.node[n]['outscore']/len(Ggraph)))
#                if float(np.log2(0.001+(Ggraph[n][n]['weight'])))>0:
                select2.append((n,G[n][n]['weight']/len(Ggraph)))
    
    ### PLOT OUTSCORE         
    selected = sorted(select, key=lambda tup: tup[1])
    tfname, outsc = zip(*selected)
    colors=[]
    place1=[]
    name1=[]
    for i in tfname:
        if i=="ZNF281" or i=="NANOG" or i=="POU5F1" or i=="SOX2" or i=="AR" or i=="ESR1" or i=="SP1": #i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="KLF7" or i=="ISL1" or i=="PHOX2B": 
            colors.append('r')
            place1.append(tfname.index(i))
            name1.append(i)
        else:
            colors.append('b')
    indices = np.arange(len(selected))
    fig = plt.figure()
    plt.hist(outsc)
    fig.savefig(imgfolder+"Distribution_Outscore_"+namefile+"_"+nameppi+".png")
    
    fig=plt.figure()
    plt.bar(indices, outsc, color=colors)
    plt.xticks(place1, name1, rotation='vertical')
    plt.tight_layout()
    fig.savefig(imgfolder+"Outscore_"+namefile+"_"+nameppi+".png")
     
    ### PLOT AUTOLINK SCORE  
    selected2 = sorted(select2, key=lambda tup: tup[1])
    tfname2, autow = zip(*selected2)
    colors2=[]
    place2=[]
    name2=[]
    for i in tfname2:
        if i=="ZNF281" or i=="NANOG" or i=="POU5F1" or i=="SOX2" or i=="AR" or i=="ESR1" or i=="SP1": #i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="KLF7" or i=="ISL1" or i=="PHOX2B": 
            colors2.append('r')
            place2.append(tfname2.index(i))
            name2.append(i)
        else:
            colors2.append('b')
    indices2 = np.arange(len(selected2))
    fig=plt.figure()
    plt.hist(autow)
    fig.savefig(imgfolder+"Distribution_Autolink_"+namefile+"_"+nameppi+".png")
    
    fig=plt.figure()
    plt.bar(indices2, autow, color=colors2)
    plt.xticks(place2, name2, rotation='vertical')
    plt.tight_layout()
    fig.savefig(imgfolder+"Autolink_"+namefile+"_"+nameppi+".png")
    
    ### PLOT INTERCONNECTIVITY SCORE
    links2 = sorted(links, key=lambda tup: tup[1])
    links=[]
    lname, allinks = zip(*links2)
    colors=[]
    place=[]
    name=[]
    for i,j in lname:
            if (i=="NANOG" and (j=="POU5F1" or j=="SOX2")) or (i=="POU5F1" and j=="SOX2"):  #i=="GATA3" or i=="HAND1" or i=="HAND2" or i=="TWIST1" or i=="MYCN" or i=="PHOX2A" or i=="PHOX2B":
                colors.append('r')
                place.append(lname.index((i,j)))
                name.append((i,j))
            else:
                colors.append('b')
    indices = np.arange(len(links2))
    fig=plt.figure()
    plt.hist(allinks)
    fig.savefig(imgfolder+"Distribution_Interconnectivity_"+namefile+"_"+nameppi+".png")
    
    fig=plt.figure()
    plt.bar(indices, allinks)
    plt.xticks(place, name, rotation='vertical')
    plt.tight_layout()
    fig.savefig(imgfolder+"Interconnectivity_"+namefile+"_"+nameppi+".png")
    
    
    score={}
    clique1=[]
    clique2=[]
    clique3=[]
    maxi1=0
    maxi2=0
    maxi3=0
    i=0
    
    print(time.ctime())
    print("finding best cliques... 6/8")
    
    ### For normalization, could be useful for subnetworks
    #Gclique=nx.find_cliques(Gc)
    #listS=[]
    #listA=[]
    #listO=[]
    #
    #for n in Gclique:
    #    i+=1
    #    listS.append((sum(G[m1][m2]['weight'] for m1 in n for m2 in n)/(len(n)*(len(n)-1))))
    #    listA.append((sum(G[m1][m1]['weight'] for m1 in n))/(len(n)))
    #    listO.append((sum(G.node[m]['outscore'] for m in n))/(len(n)))
    #print("nombre de cliques :", i)
    #
    #maxS=max(listS)
    #maxA=max(listA)
    #maxO=max(listO)
    
    Gclique=nx.find_cliques(Gc)
    ### Finding the 3 best cliques
    for n in Gclique:
        if len(n)>1:
            tmps=float(np.log2(0.001+(((sum(G[m1][m2]['weight'] for m1 in n for m2 in n))/(len(n)*(len(n)-1))))))  #mean of all edges scores
        else:
            tmps=float(np.log2(0.001+(((sum(G[m1][m2]['weight'] for m1 in n for m2 in n))/(len(n))))))
        tmpauto=((sum(G[m1][m1]['weight'] for m1 in n))/(len(n)))  #mean of Autolink scores
        tmpout=((sum(G.node[m]['outscore'] for m in n))/len(n))    #mean of Outscores
        
    #    Cscore = math.sqrt((1-tmpout)**2+(1-tmpauto)**2+(1-tmps)**2)
        Cscore=tmps*tmpauto*tmpout 
        
        ### Find a way to combine those 3 scores into a single clique score (Cscore)
        ### Probably could be improved
        
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
    
    bestclique=[clique1,clique2,clique3]
    listclique=[]
           
    listclique=bestclique  #remove this line if using subnetworks
    ###  SUBNETWORKS (not working)
    #    for n in bestclique:
    #        tmps=((sum(G[m1][m2]['weight'] for m1 in n for m2 in n))/(len(n)*(len(n)-1)))
    #        tmpauto=((sum(G[m1][m1]['weight'] for m1 in n))/(len(n)))
    #        tmpout=((sum(G.node[m]['outscore'] for m in n))/len(n))
    #        
    #    #    Cscore = math.sqrt((1-tmpout)**2+(1-tmpauto)**2+(1-tmps)**2)
    #        Cscore=tmps+tmpauto+tmpout  
    #        
    #        startS=max(G[m1][m2]['weight'] for m1 in n for m2 in n)
    #        startO=max(G.node[m]['outscore'] for m in n)
    #        startA=max(G[m1][m1]['weight'] for m1 in n)
    #        
    #        sublist=[]
    #        for (a,b) in G.edges(n):
    #            if G[a][b]['weight']==startS:
    #                sub1=a
    #                sub2=b
    #                sublist.append(sub1)
    #                sublist.append(sub2)
    #        Subout=((G.node[sub1]['outscore']+G.node[sub2]['outscore'])/2)/startO
    #        Subauto=((G[sub1][sub1]['weight']+G[sub2][sub2]['weight'])/2)/startA
    #        Subs=((G[sub1][sub2]['weight']+G[sub2][sub1]['weight'])/2)/startS
    #        
    #        Subscore = math.sqrt((1-Subout)**2+(1-Subauto)**2+(1-Subs)**2)
    #    #    Subscore=Subout+Subauto+Subs
    #        
    #        
    #        best=0
    #        cand=1
    #        while cand!=0:
    #            cand=0
    #            Subout=((sum(G.node[m]['outscore'] for m in sublist))/len(sublist))/startO
    #            Subauto=((sum(G[m1][m1]['weight'] for m1 in sublist))/len(sublist))/startA
    #            Subs=((sum(G[m1][m2]['weight'] for m1 in sublist for m2 in sublist))/(len(sublist)*(len(sublist)-1)))/startS
    #            current = math.sqrt((1-Subout)**2+(1-Subauto)**2+(1-Subs)**2)
    #    #        current = Subout+Subauto+Subs
    #            print("first",current)
    #            for tf in n:
    #                if tf not in sublist:
    #                    sublist.append(tf)
    #                    
    #                    Subout=((sum(G.node[m]['outscore'] for m in sublist))/len(sublist))/startO
    #                    Subauto=((sum(G[m1][m1]['weight'] for m1 in sublist))/len(sublist))/startA
    #                    Subs=((sum(G[m1][m2]['weight'] for m1 in sublist for m2 in sublist))/(len(sublist)*(len(sublist)-1)))/startS
    #                    print(math.sqrt((1-Subout)**2+(1-Subauto)**2+(1-Subs)**2))
    #                    if(math.sqrt((1-Subout)**2+(1-Subauto)**2+(1-Subs)**2))<current:
    #    #                if(Subout+Subauto+Subs)>current:
    #                        current=math.sqrt((1-Subout)**2+(1-Subauto)**2+(1-Subs)**2)
    #    #                    current = Subout+Subauto+Subs
    #                        best=current
    #                        cand=tf
    #                    sublist.remove(tf)
    #            if cand!=0:
    #                sublist.append(cand)
    #        Subscore=best
    #        
    #        listclique.append(sublist)                
    #    if Subscore>Cscore:
    #        listclique.append(sublist)
    #    else:
    #        listclique.append(n)
    ### SUBNETWORKS   
    
    
    
    print(time.ctime())
    print("finding most regulated genes... 7/8")    
    reg={}
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
                        
    
    ### Directed Networks corresponding to the 3 best cliques
    ### Can be exported to be visualized by gephi
    Gcl1=nx.DiGraph()
    for n in clique1:
        Gcl1.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique1:
            Gcl1.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            weight=G[n][m]['weight'])
    
    Gcl2=nx.DiGraph()
    for n in clique2:
        Gcl2.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique2:
            Gcl2.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            weight=G[n][m]['weight'])
            
    Gcl3=nx.DiGraph()
    for n in clique3:
        Gcl3.add_node(n, inscore=G.node[n]['inscore'], outscore=G.node[n]['outscore'], IsTF='core')
        G.node[n]['IsTF']='core'
        for m in clique3:
            Gcl3.add_edge(n, m, site=G[n][m]['site'], peaks=len(G[n][m]['peaks']), 
                            weight=G[n][m]['weight']) 
       
    
    ### Correlation to expression data (not working, not finished)
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
    #    plt.bar([co1,co2,co3]) 
    #    plt.show(block=False)     
    
      
    print(time.ctime())     
    print("writing files... 8/8")
    
    Gsub_path=filep.replace(".summary.bed",".gene_list")
    for n in reg:
        with open(Gsub_path+"."+str(n)+".txt","w") as file:
           for i in reg[n]:
                file.write(i+"\t"+str(reg[n][i]['score'])+"\t"+G.node[i]['IsTF']+"\t"+','.join(reg[n][i]['liste'])+"\n")
    
    nx.write_gexf(Gcl1,Gsub_path+".1.gexf")
    #nx.write_adjlist(Gcl1,Gsub_path+".1.adjlist")
    nx.write_gexf(Gcl2,Gsub_path+".2.gexf")
    #nx.write_adjlist(Gcl2,Gsub_path+".2.adjlist")
    nx.write_gexf(Gcl3,Gsub_path+".3.gexf")
    #nx.write_adjlist(Gcl3,Gsub_path+".3.adjlist")
    
    print(time.ctime())
    print("done !")

if __name__ == '__main__':
    main()