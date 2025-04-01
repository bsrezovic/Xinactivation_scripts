#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)  
library(data.table)
pathToFiles <- "/common/bsrezovic/scRNA/Wxx_2/"
output_path <- "/common/bsrezovic/scRNA/W99_2/tableWxx.csv"
SearchArea <- GRanges(seqnames = c("chrX"), ranges = IRanges(1, 156040895))
reads <- impBamGAL(pathToFiles,SearchArea,verbose = F) 
heterozygotePositions <- scanForHeterozygotes(reads, verbose =F)
countList <- getAlleleCounts(BamList = reads, GRvariants =  heterozygotePositions, verbose = F) 
a.simple <- ASEsetFromCountList(heterozygotePositions, countList)
ref(a.simple) <- randomRef(a.simple)  
genotype(a.simple) <- inferGenotypes(a.simple)
alt(a.simple) <- inferAltAllele(a.simple)
range <- granges(a.simple)
cellnames <- paste0(rep("cell_",ncol(a.simple)), as.character(c(1:ncol(a.simple))))
countnames <- paste0(rep("snpcounts_",ncol(a.simple)), as.character(c(1:ncol(a.simple))))
table <- data.table(genes = names(range),  pos = start(range), chr = as.character(seqnames(range)), ref = range$ref, alt = range$alt)
frequencies <- fraction(a.simple)
colnames(frequencies) <- cellnames
counts <- countsPerSnp(a.simple)
colnames(counts) <- countnames
tot_table <- cbind(table ,frequencies, counts)
write.csv(tot_table,output_path, row.names = T)

  
  
  