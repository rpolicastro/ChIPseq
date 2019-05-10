#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

options <- matrix(c(
	"download", "d", 1, "character", "download files from SRA? [TRUE,FALSE]",
	"samplesheet", "i", 1, "character", "required sample sheet",
	"seqdir", "o", 1, "character", "directory to download fastq files to"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## grabbing files
## --------------

if (opt$download == "TRUE") {
	# getting sample names
	samples <- read.delim(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	samples <- c(samples$R1, samples$R2, samples$R1_control, samples$R2_control)
	samples <- samples[!is.na(samples)] %>% unique
	samples <- substr(samples, 1, nchar(samples)-6)

	# downloading files
	setwd(opt$seqdir)
	for (sample in samples) {
		command <- paste("fastq-dump", sample)
		system(command)
	}
}
