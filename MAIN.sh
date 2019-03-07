#!/bin/bash

cd $PBS_O_WORKDIR
source SETTINGS

#########################################
## ChIP-seq Analysis of Multiple Samples
#########################################

### loading conda environment

source activate chipseq

### fastqc of reads

# creating directory to output fastqc results
mkdir -p ${BASEDIR}/results/fastqc

# saving fastq file names to variable
FASTQ_FILES=$(find ${SEQDIR} -name "*fastq")

# running fastqc
fastqc \
-t $CORES \
-o ${BASEDIR}/results/fastqc \
$FASTQ_FILES

### generating bowtie2 index

# create directory for genomic index
mkdir -p ${BASEDIR}/genome

# create bowtie2 index
bowtie2-build \
-f --threads $CORES \
$GENOME_FASTA \
${BASEDIR}/genome/hg38

### aligning reads to genomic index

# create directory to output alignments
mkdir -p ${BASEDIR}/results/aligned

# sample sheet aware read alignment
python ${BASEDIR}/bin/alignReads.py \
--outDir ${BASEDIR} \
--threads $CORES \
--sampleSheet $SAMPLE_SHEET

# converting to coordinate sorted bam with index
for SAM in ${BASEDIR}/results/aligned/*sam; do
	samtools sort \
	-O BAM -@ $CORES \
	-o ${BASEDIR}/results/aligned/$(basename $SAM .sam).bam \
	$SAM
done
for BAM in ${BASEDIR}/results/aligned/*bam; do samtools index $BAM; done

### calling peaks with macs2

# creating directory to output peaks
mkdir -p ${BASEDIR}/results/peaks

# calling peaks
python ${BASEDIR}/bin/callPeaks.py \
--outDir ${BASEDIR} \
--threads $CORES \
--sampleSheet $SAMPLE_SHEET
