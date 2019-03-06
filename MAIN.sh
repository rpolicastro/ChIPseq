#!/bin/bash

cd $PBS_O_WORKDIR
source SETTINGS.sh

#########################################
## ChIP-seq Analysis of Multiple Samples
#########################################

### loading conda environment

source activate chipseq

### fastqc of reads

# creating directory to output fastqc results
mkdir -p ${BASEDIR}/results/fastqc

# saving fastq files variable
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

# align reads

