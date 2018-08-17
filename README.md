# MasterTF

1. You need to have in the same directory pipeline.sh, intersect.py, enhancers.py, enhancers_left.py, filter_tads.py, merge_TADS.py, TFBS_finder.py, crossref.py, networks.py, sarus-01Mar2018.jar, filter_cream.py, creamtest.R. 
2. In the same directory you should have another directory "data", where there will be .narrowPeak and .bw files, gencode.v19.annotation.gff3, hg19_enhancer_tss_associations_FANTOM5data.bed, allTADS.bed, hg19.fa, 9606.protein.links.full.v10.5.filtered.txt, HOCOMOCOv11_full_annotation_HUMAN_mono.tsv, HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt and "pwm" directory, where you have all TF pwm files.
3. You need to have bedtools2, python3 on your computer.
4. python3 modules to install : networkX, glob, matplotlib, tqdm, numpy, pyBigWig

To run :
./pipelineCREAM.sh -i $file.narrowPeak -o outputDir -m imgDir -g gencode.v19.annotation.gff3 -t $threshold -f hg19_enhancer_tss_associations_FANTOM5data.bed -c allTADS.bed -a hg19.fa -h HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt -d HOCOMOCOv11_core_pwm_HUMAN_mono_pwm/ -b HOCOMOCOv11_full_annotation_HUMAN_mono.tsv -n True/False -p ProteinLinksFile.txt

$threshold can be calculated via a specific HMCan script

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
The pipeline begins with the gencode.v19.annotation.gff3 file. It uses it to get only transcripts, find all the TSS, and output them in the gencode.v19.TSS.txt file.
The pipeline begins with the .narrowPeak and .wig.bw files.



Optionnal :
Get a Protein-Protein Interaction file from STRING, and have it in the same format as 9606.protein.links.full.v10.5.filtered.txt.
