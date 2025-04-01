#!/usr/bin/env Rscript

library(Rsubread)
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
	  stop("At least two arguments must be supplied (folder name with bam files and an output file name for the matrix of counts)", call.=FALSE)
}
setwd(paste("./",args[1],sep=""))
results <- featureCounts(list.files(), strandSpecific = 0, annot.inbuilt = "hg38", verbose = F, nthreads = 4)

write.csv(results$counts, args[2])












