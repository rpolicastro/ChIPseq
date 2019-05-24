#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

options <- matrix(c(
	"outdir", "o", 1, "character", "output directory",
	"threads", "t", 1, "integer", "number of CPU cores",
	"samplesheet", "s", 1, "character", "required sample sheet",
	"genomesize", "g", 1, "integer", "effective genome size"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## functions
## ---------

call.peaks <- function(row) {
	command <- paste(
		"macs2 callpeak",
		"-t", file.path(opt$outdir, "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bam")),
		"-n", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"]),
		"--outdir", file.path(opt$outdir, "peaks"),
		"-g", opt$genomesize
	)
	if (!is.na(row["R2"]) & row["R2"] != "") {
		command <- paste(command, "-f BAMPE")
	} else {
		command <- paste(command, "-f BAM")
	}

	if (!is.na(row["control_ID"]) & row["control_ID"] != "") {
		command <- paste(command, "-c", file.path(opt$outdir, "aligned", paste0(row["control_ID"], ".bam")))
	}
	
	system(command)
}

## call peaks
## ----------

sample.sheet <- read.table(opt$samplesheet, header=TRUE, sep="\t", stringsAsFactors=FALSE)
apply(sample.sheet, 1, call.peaks)
