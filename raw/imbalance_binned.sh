#!/bin/bash
echo "Processing folder: $1";
j=0
for i in {1..148238850..7802000}
  do
     j=$((j+1)) 
     end=$(($i + 7802000)) 
     outfile="/common/bsrezovic/Xinactivation/validation/individual1_bams/outfiles/table_ind1_$j.csv"
     logfile="/common/bsrezovic/Xinactivation/validation/individual1_bams/logfolder/log_ind1_$j.txt"
     taskset -c 0,1,2,3,4,5,6,7 Rscript new_imbalance.R $outfile $i $end &> $logfile 
 done

