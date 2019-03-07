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
)

opt_parser  <-  OptionParser(option_list=option_list)
opt         <-  parse_args(opt_parser)

## peak annotation
## ---------------
