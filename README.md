# MasterTF

1. You need to have in the same directory pipeline.sh, intersect.py, enhancers.py, enhancers_left.py, filter_tads.py, merge_TADS.py, TFBS_finder.py, crossref.py, networks.py, sarus-01Mar2018.jar, filter_cream.py, creamtest.R. 
2. In the same directory you should have another directory "data", where will be .narrowPeak files, gencode.v19.annotation.gff3, hg19_enhancer_tss_associations_FANTOM5data.bed, allTADS.bed, hg19.fa, 9606.protein.links.full.v10.5.filtered.txt, HOCOMOCOv11_full_annotation_HUMAN_mono.tsv, HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt and "pwm" directory, where you have all TF pwm files.
3. You need to have bedtools2, python3 on your computer.
4. python3 modules to install : networkX, glob, matplotlib, tqdm, numpy, pyBigWig

To run :
./pipelineCREAM.sh -i $file.narrowPeak -o outputDir -m imgDir -g gencode.v19.annotation.gff3 -t $threshold -f hg19_enhancer_tss_associations_FANTOM5data.bed -c allTADS.bed -a hg19.fa -h HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt -d HOCOMOCOv11_core_pwm_HUMAN_mono_pwm/ -b HOCOMOCOv11_full_annotation_HUMAN_mono.tsv -n True/False -p ProteinLinksFile.txt
