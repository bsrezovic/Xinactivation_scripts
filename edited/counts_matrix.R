#!/usr/bin/env Rscript

library(Rsubread) #for feature counting, count matrix creation


args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
	  stop("At least two arguments must be supplied (folder name with bam files and an output file name for the matrix of counts)", call.=FALSE)
}
#set the directory to argument; meaning the output will go there too!
# the files must be either in the SAM or BAM folder, each representing a cell of a sample
# the files should be mapped to hg38 (or change the command bellow)
setwd(paste("./",args[1],sep=""))


# set the thread count to higher if available
results <- featureCounts(list.files(), strandSpecific = 0, annot.inbuilt = "hg38", verbose = F, nthreads = 4)


# writ the results to a .csv file (a matrix of counts)
write.csv(results$counts, args[2])












