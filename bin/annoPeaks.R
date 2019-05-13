#!/usr/bin/env Rscript

library("getopt")
library("GenomicRanges")
library("ChIPseeker")
library("rtracklayer")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("dplyr")

## command line arguments
## ----------------------

options <- matrix(c(
	"outdir", "o", 1, "character", "output directory",
	"upstream", "u", 1, "integer", "bases upstream of TSS to consider as promoter",
	"downstream", "d", 1, "integer", "bases downstream of TSS to consider as promoter"
), byrow=TRUE, ncol=5)

opt <- getopt(options)

## peak annotation
## ---------------

# ensure the extra columns in the macs2 narrowpeaks are carried over
extraCols_narrowPeak  <-  c(signal.value = "numeric", p.value.negLog10 = "numeric", q.value.negLog10 = "numeric", peak = "integer")

# go through each narrowpeaks file and annotate the peaks
for (file in list.files(file.path(opt$outdir, "peaks"), pattern=".*\\.narrowPeak")) {
	# importing narrowpeaks as GRanges object
	peaks <- import(file, format="BED", extraCols=extraCols_narrowPeak)
	# ChIPseeker to anntoate peaks using hg38 TxDB and OrgDb objects
	annotated  <-  annotatePeak(
		peaks,
		tssRegion=c(-opt$upstream, opt$downstream),
		TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
		level="transcript",
		sameStrand=FALSE,
		annoDb="org.Hs.eg.db"
	)
	
	# exporting some of the handy ChIPseeker graphs
	pdf(file.path(opt$outdir, "annotated_peaks", paste0(tools::file_path_sans_ext(basename(file)), ".pdf")))
	plotAnnoPie(annotated)
	plotAnnoBar(annotated)
	plotDistToTSS(annotated)
	dev.off()

	# writing annotated peaks to tsv
	annotated  <- as.data.frame(annotated)
	write.table(
		annotated,
		file.path(opt$outdir, "annotated_peaks", paste0(tools::file_path_sans_ext(basename(file)), ".tsv")),
		sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=""
	)
}
