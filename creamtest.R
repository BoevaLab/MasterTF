library('CREAM')

args = commandArgs(trailingOnly=TRUE)

IdentifiedCOREs <- CREAM(args[1], WScutoff=1.3, MinLength = 1000, peakNumMin = 2 )
#"/home/antoine/NAS_Public/data/projects/Antoine_Networks/results/outputDir/DNAse.H1hES.ENCODE.HMCan_peaks.narrowPeak.prepcream.bed"

#InputData <- read.table(args[1], sep="\t")
#colnames(InputData) <- c("chr", "start", "end")
#MinLength <- 1000
#if(nrow(InputData) < MinLength){
#  stop(paste( "Number of functional regions is less than ", MinLength, ".", sep = "", collapse = "")) }
#peakNumMin <- 2
#WScutoff <- 1.3
#WindowVecFinal <- WindowVec(InputData, peakNumMin, WScutoff)
#IdentifiedCOREs <- ElementRecog(InputData, WindowVecFinal, 20, peakNumMin)
#COREs <- as.data.frame(IdentifiedCOREs[c(4,2,3)])

write.table(IdentifiedCOREs, file = args[2], col.names = F, row.names = F, quote=F)