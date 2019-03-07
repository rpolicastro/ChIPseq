#!/usr/bin/env Rscript

library("getopt")
library("GenomicRanges")
library("ChIPseeker")
library("rtracklayer")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

orgdb  <-  org.Hs.eg.db
txdb   <-  TxDb.Hsapiens.UCSC.hg38.knownGene

## command line arguments
## ----------------------

option.list <- list(
	make_option(c("-o", "--outdir"), actions="store", type="character", default=getwd(), help="output directory"),
	make_options(c("-u", "--upstream"), actions="store", type="character", default=1000, help="bases upstream of TSS to consider as promoter")
	make_options(c("-d", "--downstream"), actions="store", type="character", default=1000, help="bases downstream of TSS to consider as promoter")
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## peak annotation
## ---------------

extraCols_narrowPeak  <-  c(signal.value = "numeric", p.value = "numeric", q.value = "numeric", peak = "integer")

for (file in list.files(file.path(opt$outdir, "results", "aligned"), pattern=".*\\.narrowPeak")) {
	peaks      <-  import(file, format="BED", extraCols=extraCols_narrowPeak)
	annotated  <-  annotatePeak(
		peaks,
		tssRegion=c(-opt$upstream, opt$downstream),
		TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
		level="gene",
		annoDb="org.Hs.eg.db",
	)

	pdf(file.path(opt$outdir, "results", "peaks_annotated", paste0(tools::file_path_sans_ext(basename(file)), ".tsv"))
	plotAnnoPie(annotated)
	plotAnnoBar(annotated)
	plotDistToTSS(annotated)
	dev.off()

	annotated  <- as.data.frame(annotated)
	write.table(
		annotated,
		file.path(opt$outdir, "results", "peaks_annotated", paste0(tools::file_path_sans_ext(basename(file)), ".tsv")),
		sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=""
	)
