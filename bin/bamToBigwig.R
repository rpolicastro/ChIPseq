#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line options

options <- matrix(c(
	"outdir", "d", 1, "character", "output directory",
	"samplesheet", "s", 1, "character", "required sample sheet",
	"threads", "t", 1, "integer", "number of CPU cores"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## functions

bam.to.bigwig <- function(row) {
	command <- paste(
		"bamCoverage",
		"-b", file.path(opt$outdir, "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bam")),
		"-o", file.path(opt$outdir, "bigwigs", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bigwig")),
		"-of bigwig -bs 1 --normalizeUsing CPM",
		"-p", opt$threads
	)
	if (!is.na(row["R2"]) & row["R2"] != "") {
		command <- paste(command, "-e")
	} else {
		command <- paste(command, "-e 200")
	}
	system(command)
}

control.bam.to.bigwig <- function(row) {
	command <- paste(
		"bamCoverage",
		"-b", file.path(opt$outdir, "aligned", paste0(row["control_ID"], ".bam")),
		"-o", file.path(opt$outdir, "bigwigs", paste0(row["control_ID"], ".bigwig")),
		"-of bigwig -bs 1 --normalizeUsing CPM",
		"-p", opt$threads
	)
	if (!is.na(row["R2_control"]) & row["R2_control"] != "") {
		command <- paste(command, "-e")
	} else {
		command <- paste(command, "-e 200")
	}
	system(command)
}

## bigwigs to bams

sample.sheet <- read.delim(opt$samplesheet, sep="\t", header=TRUE, stringsAsFactors=FALSE)
apply(sample.sheet, 1, bam.to.bigwig)

# dealing with control samples
if ((sample.sheet %>% filter(!is.na(control_ID) & control_ID != "") %>% nrow) > 0) {
	controls <- sample.sheet %>% distinct(control_ID, .keep_all=TRUE)
	apply(controls, 1, control.bam.to.bigwig)
}
