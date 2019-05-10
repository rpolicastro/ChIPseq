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
${OUTDIR}/genome/index/hg38

# create directory to output alignments
mkdir -p ${BASEDIR}/results/aligned

# read alignment
Rscript ${BASEDIR}/bin/alignReads.R \
--outdir $BASEDIR \
--seqdir $SEQDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET

# converting to coordinate sorted bam with index
for SAM in ${BASEDIR}/results/aligned/*sam; do
	samtools sort \
	-O BAM -@ $CORES \
	-o ${BASEDIR}/results/aligned/$(basename $SAM .sam).bam \
	$SAM
done
for BAM in ${BASEDIR}/results/aligned/*bam; do samtools index $BAM; done

## peak calling and annotation
## ---------------------------

# creating directory to output peaks
mkdir -p ${BASEDIR}/results/peaks

# calling peaks
Rscript ${BASEDIR}/bin/callPeaks.R \
--outdir $BASEDIR \
--threads $CORES \
--samplesheet $SAMPLE_SHEET

# creating directory to output annotated peak files
mkdir -p ${BASEDIR}/results/annotated_peaks

# annotating peaks
Rscript ${BASEDIR}/bin/peakAnno.R \
--outdir $BASEDIR \
--upstream $UPSTREAM \
--downstream $DOWNSTREAM

## bam to bigwig
## -------------

# creating directory to output bigwigs
mkdir -p ${BASEDIR}/results/bigwigs

# converting bams to bigwigs
Rscript ${BASEDIR}/bin/bamToBigwig.R --outdir $BASEDIR

