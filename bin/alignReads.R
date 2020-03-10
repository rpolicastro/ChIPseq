#!/usr/bin/env Rscript

library("getopt")
library("dplyr")

## command line arguments
## ----------------------

options <- matrix(c(
	"outdir", "d", 1, "character", "output directory",
	"samplesheet", "s", 1, "character", "required sample sheet",
	"seqdir", "o", 1, "character", "directory with fastq files",
	"threads", "t", 1, "integer", "number of CPU cores"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## functions
## ---------

align.experimental <- function(row) {
	# getting command ready
	command <- paste(
		"bowtie2",
		"-x", file.path(opt$outdir, "genome", "index", "genome"),
		"-S", file.path(opt$outdir, "aligned", paste0(row["sample_ID"], "_", row["condition"], "_", row["replicate"], ".sam")),
		"-q --phred33 --no-unal --threads", opt$threads
	)
	# adding paired versus unpaired options
	if (!is.na(row["R2"]) & row["R2"] != "") {
		command <- paste(command,
			"-1", file.path(opt$seqdir, row["R1"]),
			"-2", file.path(opt$seqdir, row["R2"]),
			"--no-mixed --no-discordant --dovetail -I 10 -X 700"
		)
	} else {
		command <- paste(command, "-U", file.path(opt$seqdir, row["R1"]))
	}
	# submitting command
	system(command)
}

align.control <- function(row) {
	# getting command ready
	command <- paste(
		"bowtie2",
		"-x", file.path(opt$outdir, "genome", "index", "genome"),
		"-S", file.path(opt$outdir, "aligned", paste0(row["control_ID"], ".sam")),
		"-q --phred33 --no-unal --threads", opt$threads
	)
	# adding paired versus unpaired options
	if (!is.na(row["R2_control"]) & row["R2_control"] != "") {
		command <- paste(command,
			"-1", file.path(opt$seqdir, row["R1_control"]),
			"-2", file.path(opt$seqdir, row["R2_control"]),
			"--no-mixed --no-discordant --dovetail -I 10 -X 700"
		)
	} else {
		command <- paste(command, "-U", file.path(opt$seqdir, row["R1_control"]))
	}
	# submitting command
        system(command)
}

## aligning reads
## --------------

sample.sheet <- read.table(opt$samplesheet, sep="\t", header=T, stringsAsFactors=FALSE)
apply(sample.sheet, 1, align.experimental)

# aligning control samples also
if ((sample.sheet %>% filter(!is.na(control_ID)) %>% nrow) > 0) {
	unique.controls <- sample.sheet %>% distinct(control_ID, .keep_all=TRUE)
	apply(unique.controls, 1, align.control)
}
