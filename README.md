# ChIPseq
Automation of ChIP-seq Workflow

## About

This repository will use software in a conda virtual environment to process multiple ChIP-seq samples. Before proceeding make sure you have your conda environment installed, and the settings properly set in the 'settings.conf' file.

# Getting Started

## Cloning Repository

To get started, you must first clone the ChIPseq automation repository. Navigate to a directory you would like to clone the rep to and enter `git clone https://github.com/rpolicastro/ChIPseq.git`.

## Preparing Conda Environment

This workflow takes advantage of the [conda](https://conda.io/en/latest/) package manager and virtual environment. The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty.

Before creating the environment, you must first install miniconda.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update -n base -c defaults conda`.

You are now ready to create the virtual sofware environment, and download all software and dependencies.

1. Create the new environment and specify the software to include in it.
```
conda create -n chipseq-automation -y -c conda-forge -c bioconda \
fastqc bowtie2 samtools macs2 deeptools bedtools sra-tools r-tidyverse r-getopt \
bioconductor-chipseeker bioconductor-rtracklayer bioconductor-genomicranges \
bioconductor-org.hs.eg.db bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
```
2. Update the software to the latest compatible versions.
```
conda update -n chipseq-automation -y -c conda-forge -c bioconda --all
```

If you wish to use any of the software in the environment outside of the workflow you can type `conda activate chipseq-automation`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.

## Creating Sample Sheet

In order to keep track of samples, this workflow requires the creation of a sample sheet. An example sheet 'samples.txt' is provided in the main repository directory. It is important to follow exact formatting of this sheet, as the information within it is used in various stages of the workflow.

| Column | Description |
| ------ | ----------- |
| sample_ID | Short sample identifier (e.g. A001). |
| condition | Experimental condition (e.g. EWSR1_KD). |
| replicate | Sample replicate number (e.g. 1). |
| R1 | Name of R1 fastq file of experimental condition. |
| R2 | Name of R2 fastq file of experimental condition (put NA if single end). |
| control_ID | Unique identifier for input/control sample (e.g. KD_input_1). |
| R1_control | Name of R1 fastq of input/control. |
| R2_control | Name of R2 fastq of input/control (put NA if single end). |
| paired | Put 'paired' or 'unpaired' depending on the run. |

After creating the sample sheet, set the path and file name in the 'settings.conf' file.

## Specifying Run Settings

The last step is to set a few settings in the 'settings.conf' file in the main repository directory.

| Setting | Description |
| ------- | ----------- |
| BASEDIR | The output directory for the workflow results. |
| CORES | The number of CPU cores/threads. |
| SAMPLE_SHEET | The directory and name of the sample sheet (e.g. /analysis/samples.tsv). |
| SEQDIR | The directory containing the fastq files. |
| GENOME_FASTA | The directory and name of the genome assembly fasta (e.g. /analysis/genome.fasta) |
| UPSTREAM | Bases upstream of TSS for defining promoter (e.g. 1000) |
| DOWNSTREAM | Bases downstream of TSS for defining promoter (e.g. 1000) |

## Running the Workflow

After getting the conda environment ready, the sample sheet prepared, and the settings specified, you are now ready to run the workflow. Navigate to the main directory and enter 'bash main.sh'.

###### Notes for IU Folks
If you wish to submit the workflow to a compute node, you can do so by submitting it through the TORQUE resource manager `qsub -l nodes=1:ppn=8,vmem=64gb,walltime=12:00:00 main.sh`. 'ppn' specifies the threads, and 'vmem' is the virtual memory.

# Built With

This workflow would not be possible without the great software listed below.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Read quality control.
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - Short read aligner.
- [Samtools](http://www.htslib.org/) - SAM/BAM manipulation.
- [MACS](https://github.com/taoliu/MACS) - Peak caller.
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) - Bed file manipulation.
- [deepTools](https://deeptools.readthedocs.io/en/develop/) - Read normalization and heatmap generation.
- [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) - Peak annotation.
