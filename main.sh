#!/bin/bash

cd $PBS_O_WORKDIR
source settings.conf

#########################################
## ChIP-seq Analysis of Multiple Samples
#########################################

## loading conda environment
## -------------------------

source activate chipseq-automation

## download files from SRA if necessary
## ------------------------------------

# creating sequence directory
mkdir -p $SEQDIR

# downloading from SRA if necessary
Rscript ${BASEDIR}/bin/getSRA.R \
--download $DOWNLOAD \
--samplesheet $SAMPLE_SHEET \
--seqdir $SEQDIR

## fastqc of reads
## ---------------

# creating directory to output fastqc results
mkdir -p ${OUTDIR}/fastqc

# saving fastq file names to variable
FASTQ_FILES=$(find ${SEQDIR} -name "*fastq")

# running fastqc
fastqc \
-t $CORES \
-o ${OUTDIR}/fastqc \
$FASTQ_FILES

## bowtie 2 read alignment
## -----------------------

# create directory for genomic index
mkdir -p ${OUTDIR}/genome/index

# create bowtie2 index
bowtie2-build \
-f --threads $CORES \
$GENOME_FASTA \
${OUTDIR}/genome/index/genome

# create directory to output alignments
mkdir -p ${OUTDIR}/aligned

# read alignment
Rscript ${BASEDIR}/bin/alignReads.R \
--outdir $OUTDIR \
--seqdir $SEQDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET \

# converting to coordinate sorted bam with index
for SAM in ${OUTDIR}/aligned/*sam; do
	samtools sort \
	-O BAM -@ $CORES \
	-o ${OUTDIR}/aligned/$(basename $SAM .sam).bam \
	$SAM
done
for BAM in ${OUTDIR}/aligned/*bam; do samtools index $BAM; done

## peak calling and annotation
## ---------------------------

# creating directory to output peaks
mkdir -p ${OUTDIR}/peaks

# calling peaks
Rscript ${BASEDIR}/bin/callPeaks.R \
--outdir $OUTDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET \
--genomesize $GENOME_SIZE

# creating directory to output annotated peak files
mkdir -p ${OUTDIR}/annotated_peaks

# annotating peaks
Rscript ${BASEDIR}/bin/annoPeaks.R \
--outdir $OUTDIR \
--upstream $UPSTREAM \
--downstream $DOWNSTREAM \
--genomegtf $GENOME_GTF

## bam to bigwig
## -------------

# creating directory to output bigwigs
mkdir -p ${OUTDIR}/bigwigs

# converting bams to bigwigs
Rscript ${BASEDIR}/bin/bamToBigwig.R \
--outdir $OUTDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET
