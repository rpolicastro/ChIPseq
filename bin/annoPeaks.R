#!/usr/bin/env Rscript

library("getopt")
library("GenomicRanges")
library("ChIPseeker")
library("rtracklayer")
library("GenomicFeatures")
library("dplyr")

## command line arguments
## ----------------------

options <- matrix(c(
	"outdir", "o", 1, "character", "output directory",
	"genomegtf", "a", 1, "character", "directory and file name of genomic GTF/GFF",
	"upstream", "u", 1, "integer", "bases upstream of TSS to consider as promoter",
	"downstream", "d", 1, "integer", "bases downstream of TSS to consider as promoter"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## peak annotation
## ---------------

# ensure the extra columns in the macs2 narrowpeaks are carried over
extraCols_narrowPeak  <-  c(signal.value = "numeric", p.value.negLog10 = "numeric", q.value.negLog10 = "numeric", peak = "integer")

# load GTF/GFF file as TxDb object
txdb <- makeTxDbFromGFF(opt$genomegtf)

# go through each narrowpeaks file and annotate the peaks
for (file in list.files(file.path(opt$outdir, "peaks"), pattern=".*\\.narrowPeak")) {
	file <- file.path(opt$outdir, "peaks", file)
	# importing narrowpeaks as GRanges object
	peaks <- import(file, format="BED", extraCols=extraCols_narrowPeak)
	# ChIPseeker to anntoate peaks using hg38 TxDB and OrgDb objects
	annotated  <-  annotatePeak(
		peaks,
		tssRegion=c(-opt$upstream, opt$downstream),
		TxDb=txdb,
		level="transcript",
		sameStrand=FALSE
	)
	
	# exporting some of the handy ChIPseeker graphs
	pdf(file.path(opt$outdir, "annotated_peaks", paste0(tools::file_path_sans_ext(basename(file)), ".pdf")))
	print(plotAnnoPie(annotated))
	print(plotAnnoBar(annotated))
	print(plotDistToTSS(annotated))
	dev.off()

	# writing annotated peaks to tsv
	annotated  <- as.data.frame(annotated)
	write.table(
		annotated,
		file.path(opt$outdir, "annotated_peaks", paste0(tools::file_path_sans_ext(basename(file)), ".tsv")),
		sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=""
	)
}
