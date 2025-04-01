#!/bin/bash

echo "Processing folder: $1";

#set the iterating variable
j=0
#a for loop that goes through the coordinates of chrX (hg38) 
# the loop goes from 1 to (length(chrx)-step) by step
# size of the step determines number of bins; below example results in 19 bins
for i in {1..148238850..7802000}
  do
     j=$((j+1))  #increment for output file and logfile naming
     end=$(($i + 7802000)) 
     echo "$j running chunk from $i to $end of the x chromosome"
	 #the output file name of the .csv table that will contian the allelic ratios
	 #all the tables for a single sample are then simply concatenated
     outfile="/common/bsrezovic/Xinactivation/validation/individual1_bams/outfiles/table_ind1_$j.csv"
     logfile="/common/bsrezovic/Xinactivation/validation/individual1_bams/logfolder/log_ind1_$j.txt"
     echo "File output name: $outfile"
     echo "logfile name: $logfile"
    # run the new_imblance.R script  with the three command line arguments (outfile name, bin start, bin end)
     echo "taskset -c 0,1,2,3,4,5,6,7 Rscript new_imbalance.R $outfile $i $end &> $logfile"
     taskset -c 0,1,2,3,4,5,6,7 Rscript new_imbalance.R $outfile $i $end &> $logfile 
 done

