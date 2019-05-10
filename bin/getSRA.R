#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

option.list <- list(
	make_option("--download", actions="store", type="character", help="download files from SRA? [TRUE,FALSE]")
	make_option("--samplesheet", actions="store", type="character", help="required sample sheet"),
	make_option("--seqdir", acitons="store", type="character", help="directory with fastq files")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## grabbing files
## --------------

if (opt$download == "TRUE") {
	# getting sample names
	samples <- read.delim(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	samples <- c(samples$R1, samples$R2, samples$R1_control, samples$R2_control)
	samples <- samples[!is.na(samples)] %>% unique
	samples <- substr(samples, 1, nchar(samples)-6)

	# downloading files
	setwd(seqdir)
	for (sample in samples) {
		command <- paste("fastq-dump", sample)
		system(command)
	}
}
