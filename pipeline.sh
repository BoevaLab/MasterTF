#/usr/bin/env bash

while getopts i:o:g:t:f:c:a:h:d:b: option
do
case "${option}"
in
o) OUTDIR=${OPTARG};;
i) InputFile=${OPTARG};;
g) GencodeFile=${OPTARG};;
t) threshold=${OPTARG};;
f) FANTOM5=${OPTARG};;
c) TADS=${OPTARG};;
a) GENOMEFA=${OPTARG};;
h) HOCOMOCOthresholds=${OPTARG};;
d) PWMDIR=${OPTARG};;
b) HOCOMOCOtoTF=${OPTARG};;
esac
done

display_usage() { 
echo -e "\nUsage: pipeline.sh -i InputDir -o OutputDir -g gencodeFile -t threshold -f Fantom5PromEnhancerInteractionsFile -c fileWithTADs -a genome.fa -h HOCOMOCOthresholds -d folderWithPWMs -b HOCOMOCOfullannotations\n"
echo -e "\nExample: pipeline.sh -i ATAC.XXX.rep1_peaks.narrowPeak -o /home/user/outputDir -g gencode.v19.annotation.gff3 -t 4.99 -f hg19_enhancer_tss_associations_FANTOM5data.bed -c allTADS.bed -a hg19.fa -h HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt -d HOCOMOCOv11_core_pwm -b HOCOMOCOv11_full_annotation_HUMAN_mono.tsv\n"
echo "You can get all necessary HOCOMOCO files (PWMs and thresholds) here: "
echo "PWMs: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwm_HUMAN_mono.tar.gz "
echo "Thresholds: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt "
echo "HOCOMOCO to TF: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv " 
} 

SCRIPTPATH=`dirname $0`

if [ ! -f "$HOCOMOCOtoTF" ]
then
echo -e "\nERROR: You have to set the file with HOCOMOCO full annotation (e.g. hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv)"
display_usage
exit 1
fi

if [ ! -f "$HOCOMOCOthresholds" ]
then
echo -e "\nERROR: You have to set the file with HOCOMOCO thresholds (e.g. http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_standard_thresholds_HUMAN_mono.txt)"
display_usage
exit 1
fi

if [ ! -d "$PWMDIR" ]
then
echo -e "\nERROR: You have to set the directory with HOCOMOCO PWM (e.g. get it here: http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_pwm_HUMAN_mono.tar.gz)"
display_usage
exit 1
fi

if [ ! -f "$FANTOM5" ]
then
echo -e "\nERROR: You have to set the .bed file with promoter-enhancer interactions from FANTOM5"
display_usage
exit 1
fi

if [ ! -f "$GENOMEFA" ]
then
echo -e "\nERROR: You have to set the .fasta file with the genomic sequences"
display_usage
exit 1
fi


if [ ! -f "$TADS" ]
then
echo -e "\nERROR: You have to set the .bed file with TADs"
display_usage
exit 1
fi

if [ ! -f "$InputFile" ]
then
echo -e "\nERROR: You have to set the .bed file with ATAC-seq peaks (e.g. narrowPeak file from the HMCan output)"
display_usage
exit 1
fi

if [ -z "$threshold" ]
then
echo -e "\nERROR: You have to set the score threshold for the .bed file with ATAC-seq peaks"
display_usage
exit 1
fi

if [ ! -f "$GencodeFile" ]
then
echo -e "\nERROR: You have to set the file with gene annotation (.gff)"
display_usage
exit 1
fi

if [ -z "$OUTDIR" ]
then
echo -e "\nERROR: You have to set the output directory"
display_usage
exit 1
fi

echo "TFhubsFinder -- input parameters:"
echo "Input file: $InputFile"
echo "Threshold: $threshold"
echo "Output Directory: $OUTDIR"
echo "GFF file with gene information: $GencodeFile"
echo "BED file with promoter-enhancer interactions: $FANTOM5"
echo "BED file with TADs: $TADS"
echo "FASTA file with the genome: $GENOMEFA"
echo "File with HOCOMOCO thresholds: $HOCOMOCOthresholds"
echo "File with full HOCOMOCO annotations: $HOCOMOCOtoTF"

echo "TFhubsFinder create the output directory if it does not exist"

mkdir $OUTDIR || true

rm -rf $OUTDIR/FASTA/ || true

mkdir $OUTDIR/FASTA/ || true

rm -rf $OUTDIR/TFBS/ || true

mkdir $OUTDIR/TFBS/ || true

output_fasta=$OUTDIR/FASTA/
TFBS_output=$OUTDIR/TFBS/

mkdir $TFBS_output/promoters
mkdir $TFBS_output/FANTOM5
mkdir $TFBS_output/TADS


file=$(basename $InputFile)
prepared_output=$OUTDIR/${file}.prepared.bed
promoters_output=$OUTDIR/${file}.promoters.bed
uniq_promoters_output=$OUTDIR/${file}.uniq_promoters.bed
enhancers_output=$OUTDIR/${file}.enhancers.bed
FANTOM5_enhancers_output=$OUTDIR/${file}.FANTOM5_enhanc.bed
uniq_FANTOM5_enhanc_output=$OUTDIR/${file}.uniq_FANTOM5_enhanc.bed
enhancers_left_output=$OUTDIR/${file}.enhanc_left.bed
combined_enhanc_tads_output=$OUTDIR/${file}.comb_enhanc_tads.bed
unfiltered_tads_enhanc_output=$OUTDIR/${file}.unfilt_enhanc_tads.bed
TADS_with_repeats_output=$OUTDIR/${file}.TADS_with_repeats
TADS_enhancers_output=$OUTDIR/${file}.TADS_enhanc.bed
promoters_fasta_output=$output_fasta/${file}.promoters.fa
FANTOM5_enhanc_fasta_output=$output_fasta/${file}.FANTOM5.fa
TADS_enhanc_fasta_output=$output_fasta/${file}.TADS.fa
graph_output=$OUTDIR/${file}.gefx
cross_output=$OUTDIR/${file}.cross.bed

rm $OUTDIR/${file}.cross.bed.idx

date
echo "0. getting transcripts..."
tail -n +8 $GencodeFile | awk '
BEGIN {OFS="\t"}{
	if ($7 == "+") {print $1,$3,$4,$4+1,$9}
	else if ($7 == "-") {print $1,$3,$5-1,$5,$9}
}' | grep -wE "(transcript)" | awk 'BEGIN {OFS="\t"}{print $1,$2,$3,$4,$5}' > $OUTDIR/gencode.v19.TSS.txt 
echo "done!"

echo "determining promoter regions..."
cat $OUTDIR/gencode.v19.TSS.txt | awk 'BEGIN {OFS="\t"}{print $1,$3-750,$4+750,"TSS",$5}' > $OUTDIR/promoter_regions.txt.bed
echo "done!"     

echo "1. filtering $file with threshold $threshold..."
cat $InputFile | awk "BEGIN {OFS=\"\t\"} {if (\$5 > $threshold) {print \$1,\$2,\$3,\"peak\",\$4,\$10}}" > $prepared_output
#}

echo "2. looking for promoters..."
cat $prepared_output $OUTDIR/promoter_regions.txt.bed | sort -k1,1 -k2,2n | python3 $SCRIPTPATH/intersect.py TSS > $promoters_output
# output: chrom, start, end, peak, gene

echo "3. looking for enhancers..."
cat $promoters_output | sort -k1,1 -k2,2n -k5,5 | sort -u -k4,4 > $uniq_promoters_output
python3 $SCRIPTPATH/enhancers.py $uniq_promoters_output $prepared_output > $enhancers_output
  
# output: chrom, start, end, "peak", peak

echo "4. looking for fantom5 enhancers..."
cat $FANTOM5 | tail -n +3 | awk 'BEGIN {OFS="\t"} {print $1,$7-200,$7+200,"genes",$4}' > $OUTDIR/hg19_FANTOM5data.bed
cat $enhancers_output $OUTDIR/hg19_FANTOM5data.bed | sort -k1,1 -k2,2n | python3 $SCRIPTPATH/intersect.py genes > $FANTOM5_enhancers_output
cat $FANTOM5_enhancers_output | sort -k1,1 -k2,2n -k5,5 | sort -u -k4,4 > $uniq_FANTOM5_enhanc_output
echo "done!"
        
echo "5. looking for TADS enhancers..."
python3 $SCRIPTPATH/enhancers_left.py $FANTOM5_enhancers_output $enhancers_output > $enhancers_left_output
cat $TADS | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"TADS",$4}' > $combined_enhanc_tads_output
cat $enhancers_left_output $combined_enhanc_tads_output | sort -k1,1 -k2,2n | python3 $SCRIPTPATH/intersect.py TADS > $unfiltered_tads_enhanc_output
python3 $SCRIPTPATH/filter_tads.py $promoters_output $unfiltered_tads_enhanc_output >  $TADS_with_repeats_output
python3 $SCRIPTPATH/merge_TADS.py $TADS_with_repeats_output > $TADS_enhancers_output

echo "6. making fasta of promoters and enhancers..."
bedtools getfasta -fi $GENOMEFA -bed $uniq_promoters_output -name -fo $promoters_fasta_output
bedtools getfasta -fi $GENOMEFA -bed $uniq_FANTOM5_enhanc_output -name -fo $FANTOM5_enhanc_fasta_output
bedtools getfasta -fi $GENOMEFA -bed $TADS_enhancers_output -name -fo $TADS_enhanc_fasta_output
echo "done!"
        
tail -n +2 $HOCOMOCOthresholds > $OUTDIR/HOCOMOCO_thresholds.txt
date

echo "looking for TFBS in promoters..."
python3 $SCRIPTPATH/TFBS_finder.py $OUTDIR/HOCOMOCO_thresholds.txt $PWMDIR $promoters_fasta_output promoters $HOCOMOCOtoTF

echo "looking for TFBS in FANTOM5 enhancers..."
python3 $SCRIPTPATH/TFBS_finder.py $OUTDIR/HOCOMOCO_thresholds.txt $PWMDIR $FANTOM5_enhanc_fasta_output FANTOM5 $HOCOMOCOtoTF
        
echo "looking for TFBS in TADS enhancers..."
python3 $SCRIPTPATH/TFBS_finder.py $OUTDIR/HOCOMOCO_thresholds.txt $PWMDIR $TADS_enhanc_fasta_output TADS $HOCOMOCOtoTF
         
echo "7. building the network..."
date
python3 $SCRIPTPATH/crossref.py $promoters_output $FANTOM5_enhancers_output $TADS_enhancers_output $TFBS_output
python3 $SCRIPTPATH/networks.py $cross_output
echo "done!"

sort -k1,1 -k2n,2n -k4,4 $cross_output -o $cross_output

LC_ALL=C sort -k2,2rg $OUTDIR/${file}.gene_list.1.txt -o $OUTDIR/${file}.gene_list.1.txt
LC_ALL=C sort -k2,2rg $OUTDIR/${file}.gene_list.2.txt -o $OUTDIR/${file}.gene_list.2.txt
LC_ALL=C sort -k2,2rg $OUTDIR/${file}.gene_list.3.txt -o $OUTDIR/${file}.gene_list.3.txt

echo "The Network files $graph_output 1, 2 and 3 are complete, you can visualize them with Gephi"
echo "The regulated genes files ${file}.gene_list.txt 1, 2 and 3 are complete"
date

########
echo "8. comparing to Coltron..."
python3 $SCRIPTPATH/stats.py $file $OUTDIR
