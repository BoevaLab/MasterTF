# MasterTF

1. You need to have in the same directory pipeline.sh, intersect.py, enhancers.py, enhancers_left.py, filter_tads.py, merge_TADS.py, TFBS_finder.py, crossref.py, networks.py, sarus-01Mar2018.jar, filter_cream.py, creamtest.R. 
2. In the same directory you should have another directory "data", where there will be .narrowPeak and .bw files, gencode.v19.annotation.gff3, hg19_enhancer_tss_associations_FANTOM5data.bed, allTADS.bed, hg19.fa, 9606.protein.links.full.v10.5.filtered.txt, HOCOMOCOv11_full_annotation_HUMAN_mono.tsv, HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt and "pwm" directory, where you have all TF pwm files.
3. You need to have bedtools2, python3 on your computer.
4. You need to have the CREAM R package installed https://github.com/bhklab/CREAM.
5. python3 modules to install : networkX, glob, matplotlib, tqdm, numpy, pyBigWig

To run :
./pipelineCREAM.sh -i $file.narrowPeak -o outputDir -m imgDir -g gencode.v19.annotation.gff3 -t $threshold -f hg19_enhancer_tss_associations_FANTOM5data.bed -c allTADS.bed -a hg19.fa -h HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt -d HOCOMOCOv11_core_pwm_HUMAN_mono_pwm/ -b HOCOMOCOv11_full_annotation_HUMAN_mono.tsv -q 3 -n True/False -p ProteinLinksFile.txt

-q : number of best cliques found by the algorithm, -n : True or False, using the Protein-Protein Interaction option or not 

-t $threshold can be calculated via a specific HMCan script

The BigWig file must have the exact same name as the narrowPeak file, for example :
ATAC.H1hES.50Kcells.rep1.bwa.q20.rmdup_peaks.narrowPeak  and
ATAC.H1hES.50Kcells.rep1.bwa.q20.rmdup_peaks.narrowPeak.wig.bw

Files used :
Geo accession :

ATAC-Seq K562 cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1782764
  replicate 1: SRR2085918.sra
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1782765
  replicate 2: SRR2085919.sra

ATAC-Seq G12878 cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1155957
  replicate 1: SRR891268.sra
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1155958
  replicate 2: SRR891269.sra
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1155959
  replicate 3: SRR891270.sra

ATAC-Seq H1hES cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2264818
  replicate 1: SRR4011946.sra
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2264819
  replicate 2: SRR4011947.sra

ATAC-Seq NGP cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2113546

ATAC-Seq KELLY cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2113545

ATAC-Seq BE2C cells :

  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2113544


METHOD :

Before the beginning, the gencode.v19.annotation.gff3 file is used to get only transcripts, find all the transcripts TSS, and output them in the gencode.v19.TSS.txt file.

From this file are deduced all the promoter regions, they are the regions from -750bp to +750bp around the TSS.

The pipeline begins then with the .narrowPeak file. It is first filtered using the threshold parameter, only the 'best' peaks make it into the .prepared.bed file.

CREAM is used on this file to find COREs.

The promoter regions are used with the .prepared.bed file to find only the promoter peaks and output them in the .promoters.bed file, associated to the corresponding genes.

The rest of the peaks are considered enhancers. 
The ones that are in the hg19_enhancer_tss_associations_FANTOM5data.bed file are considered fantom5 enhancer peaks, and can be found in the .FANTOM5_enhanc.bed file, associated to the corresponding genes.
The rest of the peaks are considered TADS enhancers, and are associated to all genes in the same TAD, using the allTADS.bed file. They can be found in the .TADS_enhanc.bed file.

bedtools is used to make 3 fasta files corresponding to the promoters, fantom5 enhancers, and TAD enhancers.

SARUS is then used on these fasta files via the TFBS_finder.py script to do a motif analysis. The results are in the 3 folders (promoters, fantom5 and TADS), one file by Transcription Factor Motif.

The crossref.py script reads all these results, puts in one dictionary the peaks and the genes they are associated to, and in another dictionary the peaks and the TF motifs found inside them. These two informations are crossed in the .summary.bed file, where you can find which TF has a motif in which peak, in what position, with what probability of actually binding, and which genes are associated to this peak. If the probability of binding is too low (<0.1), the motif is disregarded.

The networks.py script reads this summary file and builds a network, in which the nodes are genes and TFs, and the edges are the regulations between them. There is an edge between one TF and one gene (or another TF) if the binding motif of the TF is found in a peak associated to the gene. At this step you can use the Protein-Protein Interaction option to find links between a TF and a gene, even if the TF doesn't bind directly the DNA but forms a complex with another TF. Scores are then calculated for each edge and each node.

A first selection is the used to find a list of the best TFs based on this scores, and then a clique analysis is used to find all the cliques of TFs within this list. The scores are plotted, and the X best cliques (-q parameter) are found. A list of the genes most regulated by the TFs of the clique is also generated.

Optional :

You can get a Protein-Protein Interaction file from STRING, and have it in the same format as 9606.protein.links.full.v10.5.filtered.txt. By setting the "-p" parameter to this file, and the "-n" parameter to True, the pipeline can now take into account non-direct TF binding to the DNA.

OUTPUT :

In the outputDir folder :

.gexf files, corresponding to the best cliques. Those files can be visualized using different softwares like Gephi https://gephi.org/.

.gene_list.txt, corresponding to the genes most regulated by the TFs of the clique.

.summary.bed file, corresponding to all motifs found and their scores. This file can be visualized using IGV https://software.broadinstitute.org/software/igv/.

In the imgDir folder (-m parameter) :

6 .png files, corresponding to the 3 scores calculated for the cliques (means of outscore, autolink score and interconnectivity), and their distributions.
