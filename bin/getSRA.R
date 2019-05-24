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
	if (!is.na(row["R2"]) & row["R2"] != "") {
		sra.id <- row["R1"] %>% substr(., 1, nchar(.)-8)
		command <- paste("fastq-dump", "--split-files", sra.id)
		system(command)
	} else {
		sra.id <- row["R1"] %>% substr(., 1, nchar(.)-6)
		command <- paste("fastq-dump", sra.id)
		system(command)
	}
}

grab.sra.controls <- function(row) {
	# download controls if they exist and split file if paired end
	if (!is.na(row["R2_control"]) & row["R2_control"] != "") {
		sra.id <- row["R1_control"] %>% substr(., 1, nchar(.)-8)
		command <- paste("fastq-dump", "--split-files", sra.id)
		system(command)
	} else {
		sra.id <- row["R1_control"] %>% substr(., 1, nchar(.)-6)
		command <- paste("fastq-dump", sra.id)
		system(command)
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
	
	# downloading controls
	controls <- samples %>% distinct(control_ID, R1_control, R2_control)
	apply(controls, 1, grab.sra.controls)
}
