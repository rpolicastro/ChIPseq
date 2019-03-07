#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

option.list <- list(
	make_option(c("-o", "--outdir"), actions="store", type="character", default=getwd(), help="output directory"),
	make_option("--threads", actions="store", type="numeric", default=1, help="number of CPU cores"),
	make_option("--samplesheet", actions="store", type="character", help="required sample sheet"),
	make_option("--seqdir", acitons="store", type="character", help="directory with fastq files")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## functions
## ---------

align.experimental <- function(row) {
	if (row["paired"] == "paired") {
		command <- paste(
			"bowtie2",
			"-x", file.path(opt$outdir, "genome", "hg38"),
			"-1", file.path(opt$seqdir, row['R1']),
			"-2", file.path(opt$seqdir, row['R2']),
			"-S", file.path(opt$outdir, "results", "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".sam")),
			"-q --phred33 --no-mixed --no-discordant --threads", opt$threads
		)
	} else {
		command <- paste(
                        "bowtie2",
                        "-x", file.path(opt$outdir, "genome", "hg38"),
                        "-U", file.path(opt$seqdir, row['R1']),
                        "-S", file.path(opt$outdir, "results", "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".sam")),
                        "-q --phred33 --threads", opt$threads
		)
	}
	system(command)
}

align.control <- function(row) {
        if (row["paired"] == "paired") {
                command <- paste(
                        "bowtie2",
                        "-x", file.path(opt$outdir, "genome", "hg38"),
                        "-1", file.path(opt$seqdir, row['R1_control']),
                        "-2", file.path(opt$seqdir, row['R2_control']),
                        "-S", file.path(opt$outdir, "results", "aligned", paste0(row["control_ID"], ".sam")),
                        "-q --phred33 --no-mixed --no-discordant --threads", opt$threads
                )
        } else {
                command <- paste(
                        "bowtie2",
                        "-x", file.path(opt$outdir, "genome", "hg38"),
                        "-U", file.path(opt$seqdir, row['R1_control']),
                        "-S", file.path(opt$outdir, "results", "aligned", paste0(row["control_ID"], ".sam")),
                        "-q --phred33 --threads", opt$threads
                )
        }
        system(command)
}

## aligning reads
## --------------

sample.sheet <- read.table(file.path(opt$outdir, opt$samplesheet), sep="\t", header=T)
apply(sample.sheet, 1, align.experimental)
unique.controls <- sample.sheet %>%
	select(control_ID, R1_control, R2_control, paired) %>%
	distinct(control_ID, .keep_all=TRUE)
apply(unique.controls, 1, align.control)
