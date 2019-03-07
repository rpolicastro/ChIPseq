#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

option.list <- list(
	make_option(c("-o", "--outdir"), actions="store", type="character", default=getwd(), help="output directory"),
	make_option("--threads", actions="store", type="numeric", default=1, help="number of CPU cores"),
	make_option("--samplesheet", actions="store", type="character", help="required sample sheet"),
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## functions
## ---------

call.peaks <- function(row) {
	if (row["paired"] == "paired") {
		command <- paste(
			"macs2 callpeak",
			"-t", file.path(opt$outdir, "results", "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bam")),
			"-c", file.path(opt$outdir, "results", "aligned", paste0(row["control_ID"], ".bam")),
			"-n", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"]),
			"-o", file.path(opt$outdir, "results", "peaks"),
			"-f BAMPE -g 3e9"
		)
	} else {
		command <- paste(
			"macs2 callpeak",
			"-t", file.path(opt$outdir, "results", "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".bam")),
			"-c", file.path(opt$outdir, "results", "aligned", paste0(row["control_ID"], ".bam")),
			"-n", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"]),
			"-o", file.path(opt$outdir, "results", "peaks"),
			"-f BAM -g 3e9"
	}
	system(command)
}

## call peaks
## ----------

sample.sheet <- read.table(opt$samplesheet, header=T, sep="\t")
apply(sample.sheet, 1, call.peaks())
