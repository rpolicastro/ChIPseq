#!/usr/bin/env Rscript

library("getopt")
library("tidyverse")

## command line arguments
## ----------------------

options <- matrix(c(
	"download", "d", 1, "character", "download files from SRA? [TRUE,FALSE]",
	"samplesheet", "i", 1, "character", "required sample sheet",
	"seqdir", "o", 1, "character", "directory to download fastq files to"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## functions
## ---------

grab.sra <- function(row) {
	# download R1 and split file if paired end
	if (!is.na(row["R2"])) {
		sra.id <- row["R1"] %>% substr(., 1, nchar(.)-8)
		command <- paste("fastq-dump", "--split-files", sra.id)
		system(command)
	} else {
		sra.id <- row["R1"] %>% substr(., 1, nchar(.)-6)
		command <- paste("fastq-dump", sra.id)
		system(command)
	}

	# download control
	downloaded.controls <- c()
	if (!is.na(row["control_ID"]) & !(row["control_ID"] %in% downloaded.controls)) {
		if (!is.na(row["R2_control"])) {
			sra.id <- row["R1_control"] %>% substr(., 1, nchar(.)-8)
			command <- paste("fastq-dump", "--split-files", sra.id)
			system(command)
			downloaded.controls <- c(downloaded.controls, row["control_ID"])
			
		} else {
			sra.id <- row["R1_control"] %>% substr(., 1, nchar(.)-6)
			command <- paste("fastq-dump", sra.id)
			downloaded.controls <- c(downloaded.controls, row["control_ID"])
		}
	}
}

## grabbing files
## --------------

if (opt$download == "TRUE") {
	# getting sample names
	samples <- read.delim(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% as_tibble
	
	# downloading files
	setwd(opt$seqdir)
	apply(samples, 1, grab.sra)
}
