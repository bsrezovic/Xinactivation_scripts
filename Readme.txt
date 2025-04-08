Folder containing scripts relevant to the X inactivation pulmonary hypertension paper
availible at:

The "raw" folder contains scripts without comments, "edited" has files with explanatory comments.

Procedures by files in order:

1. For allelic ratios:

	align.R 
		Input is a folder of raw read files in fastq.gz format
		This scrip aligns them to hg38
		Output is a folder of .bam files (and various metadata files)
	
	imbalance_binned.sh
		Bash script that runs sequential instances of
		new_imbalance.R on bins of chrX for memory saving purposes
		inputs are a folder of .bam files and the new_imbalance.R script
		
	new_imbalance.R
		uses the R modified R package AllelicImbalance
		input is an output file location, a set of .bam files and genomic coordinates
		output is a list of heterozygous sites and allelic ratios within those coordinates in 
		 table format
	Imbalance_analysis.Rmd
		Inputs are the tables from new_imbalance.R
 			Additional filtering done using scRNA data explained below under "Analysis_scRNA"
		outputs are graphs and lists of escapee genes

2. For single cell RNA analysis
	count_matrix.R
		Input is a folder of .bam files
		Outputs a matrix of gene counts in .csv format
	Analysis_scRNA.Rmd
		Inputs are the matrices of counts
		Creates SingleCellExperiment objects
                        An unfiltered verion of this object was sent to Ayu Hutami Syarif
			 for filtering and the results used in Imbalance_analysis.Rmd
		Filtering for mitochondrial genes etc., other quality controls
		UMAPS and other graphs
