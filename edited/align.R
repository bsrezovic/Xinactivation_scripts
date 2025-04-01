#!/usr/bin/env Rscript
library(Rsubread)
args = commandArgs(trailingOnly=TRUE)

# short if statement in case an argument is mistakenly not supplied correctly
if (length(args)==0) {
  stop("At least two arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.bam"
}

# Index is a indexed copy of hg38 reference genome
align(index = args[1], readfile1 = args[2], type = "rna", input_format = "gzFASTQ", output_format = "BAM", ,output_file = args[3] )
