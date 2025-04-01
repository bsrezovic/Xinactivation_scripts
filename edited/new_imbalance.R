#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)  # add command line arguments (needed if splitting the X chromosome into bins due to RAM requirements)
library(AllelicImbalance) # note that this one causes a lot of packages to load as requirements, i've listed them in previous emails
library(data.table)

# This script works on a folder filled with .bam and .bai files for all cells in a
# sample, and outputs a .csv file with allelic ratios for heterozygous sites 


# Folder with the .bam and .bai files
pathToFiles <- "/common/bsrezovic/scRNA/Wxx_2/"

# Output filepath for the final table
output_path <- "/common/bsrezovic/scRNA/W99_2/tableWxx.csv"

# The target area in the genome for which to scan for heterozygous sites
#  in all the cells ( the entire chromosome X - we filter par regions later).
# NOTE: for larger samples (more cells in a sample) or computer with less RAM this might be a problem,
# as AllelicImbalance requires all the .bam files for a single sample to be read into RAM at the same time in order
# to call the allelic ratios for the sites. A workaround is to modify this code to go through the coordinates for 
# chromosome X in smaller increments instead of all at once, and then just merge the resulting tables afterwards
#         simply add the desired coordinates as command line arguments and plug them into the 
#           command bellow while looping through the chunks of chrX in a bash script
#    example of such a  command: 
#    SearchArea <- GRanges(seqnames = c("chrX"), ranges = IRanges(as.numeric(args[2]), as.numeric(args[3])))
#     where arg[2] is the start of the bin and arg[3] the end (arg[1] is outfile name)
#    the script file imbalance_bined.sh shows how to generate the bins

# if you have enough RAM then the entire chromosome can be searched for heteozygous sites at onece like this:
SearchArea <- GRanges(seqnames = c("chrX"), ranges = IRanges(1, 156040895))


# Imports the specified genomic region from all the .bam files in the target folder
# and forms a 'GAlignmentsList' type object
reads <- impBamGAL(pathToFiles,SearchArea,verbose = F) # about 30 sec per .bam file

# Identifies the SNP positions 
heterozygotePositions <- scanForHeterozygotes(reads, verbose =F)# about 170 secs per .bam file

# Returns the allele counts at the discovered heterozygote positions
countList <- getAlleleCounts(BamList = reads, GRvariants =  heterozygotePositions, verbose = F) 

# Turn the countlist into a ASEset object (so it can store genotype information)
a.simple <- ASEsetFromCountList(heterozygotePositions, countList)


# assing one of the two top expressed alleles as the "reference allele"
# this function is poorly named as it doesent do this randomly, but goes by 
# allele count information, basically choosing the one that is more expressed.
# In our case this should always be the one expressed from the active X.
ref(a.simple) <- randomRef(a.simple)  

# modify the ASEset object to hold full genotype information
genotype(a.simple) <- inferGenotypes(a.simple)
# infer alt allele (does the same thing as randomRef but ignores the allele chosen as reference)
alt(a.simple) <- inferAltAllele(a.simple)

# save a GenomicRange object  of this for later use
range <- granges(a.simple)


#fixing the cell names to be of the form cell_1, cell_2 etc.
cellnames <- paste0(rep("cell_",ncol(a.simple)), as.character(c(1:ncol(a.simple))))
#fixing the cell names here now instead
countnames <- paste0(rep("snpcounts_",ncol(a.simple)), as.character(c(1:ncol(a.simple))))

# creating a table with the genomic position, reference and alt allele
table <- data.table(genes = names(range),  pos = start(range), chr = as.character(seqnames(range)), ref = range$ref, alt = range$alt)

# Calculate the allelic frequencies and store the snp count data for later filtering
frequencies <- fraction(a.simple)
colnames(frequencies) <- cellnames
counts <- countsPerSnp(a.simple)
colnames(counts) <- countnames

# final table containing genomic position, reference and alt allele, frequency and counts per snp data
tot_table <- cbind(table ,frequencies, counts)

# save it in .csvform
write.csv(tot_table,output_path, row.names = T)

  
  
  