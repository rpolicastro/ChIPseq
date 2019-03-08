#!/usr/bin/env Rscript

library("getopt")

## command line options

option.list <- list(
	make_option(c("-o", "--outdir"), actions="store", type="character", default=getwd(), help="output directory"),
	make_option(c("-s", "--samplesheet"), actions="store", type="character", help="directory and name of sample sheet"),
	make_option(c("-t", "--threads"), actions="store", type="numeric", default=1, help="number of CPU cores/threads")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## functions

bam.to.bigwig <- function(row) {
	command <- paste(
		"bamCoverage",
		"-b", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bam"),
		"-o", file.path(opt$outdir, "results", "bigwigs", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bigwig"))
		"-of bigwig -bs 1 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC",
		"-p", opt$threads
	)
	if (row["paired"] == "paired") {
		command <- paste(command, "-e")
	} else {
		command <- paste(command, "-e 200")
	}
	system(command)
}

## bigwigs to bams

sample.sheet <- read.delim(opt$samplesheet, sep="\t", header=TRUE)
apply(sample.sheet, 1, bam.to.bigwig())
