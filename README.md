# ChIPseq
Automation of ChIP-seq Workflow

## About

This repository will use software in a conda virtual environment to process multiple ChIP-seq samples. Before preceding make sure you have your conda environment installed, and the settings properly set in the 'settings.conf' file.

# Getting Started

## Cloning Repo

To get started, you must first clone the ChIPseq automation repository. Navigate to a directory you would like to clone the rep to and enter `git clone https://github.com/rpolicastro/ChIPseq.git`.

## Preparing Conda Environment

This workflow takes advantage of the [conda](https://conda.io/en/latest/) package manager and virtual environment. The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty. Furthermore, the provided 'environment.yml' file will reproduce the software environment used when developing the workflow. This ensures prolonged compatabilty and reproducibility.

Before creating the environment, you must first install miniconda on your computer.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update conda`.

You are now read to create the virtual sofware environment, and download all software and dependencies. If you would like to recreate the environment used when writing the original workflow, navigate to the main repository directory and enter `conda create -f environment.yml`. If you would like to create your own environment with the latest software version, follow the steps below.

1. Create the new environment and specify the software to include in it `conda create -n chipseq-automation -y -c conda-forge -c bioconda fastqc bowtie2 samtools macs2 pandas`.
2. Update the software to the latest compatible versions `conda update -n chipseq-automation -y -c conda-forge -c bioconda --all`.

If you wish to use any of the software in the environment outside of the workflow you can type `conda activate chipseq-automation`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.

## Assembling Sample Sheet

In order to keep track of samples, this workflow requires the creation of a sample sheet. An example sheet 'samples.txt' is provided in the main directory. It is important to follow exact formatting of this sheet, as the information within it is used in various stages of the workflow.

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

The last step is to set a few settings in the 'settings.conf' file in the main directory.

| Setting | Description |
| ------- | ----------- |
| BASEDIR | The output directory for the workflow results. |
| CORES | The number of CPU cores/threads. |
| SAMPLE_SHEET | The directory and name of the sample sheet (e.g. /analysis/samples.tsv). |
| SEQDIR | The directory containing the fastq files. |
| GENOME_FASTA | The directory and name of the genome assembly fasta (e.g. /analysis/genome.fasta) |
