# ChIPseq
Automation of ChIP-seq Workflow

## About

This repository will use software in a conda virtual environment to process multiple ChIP-seq samples. Before proceeding make sure you have your conda environment installed, and the settings properly set in the 'settings.conf' file.

# Getting Started

## Cloning Repository

To get started, you must first clone the ChIPseq automation repository. Navigate to a directory you would like to clone the repository to and enter `git clone https://github.com/zentnerlab/ChIPseq.git`.

## Preparing Conda Environment

This workflow takes advantage of the [conda](https://conda.io/en/latest/) package manager and virtual environment. The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty. Additionally, the provided 'environments.yml' file can be used to install the same major software versions as used to develop the script, ensuring prolonged compatability and reproducibility.

Before creating the environment, you must first install miniconda.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update -n base -c defaults conda`.

You are now ready to create the virtual sofware environment, and download all software and dependencies. It is recommended to reproduce the environment used when creating the script, but instructions on installing the latest software are provided as an alternative.

#### Reproducing the Development Environment (Recommended)

To install the major software versions used when developing this script, navigate to the 'DOCS' directory, and use the provided 'environments.yml' file to create your conda environment.
```
conda env create -f environment.yml
```
For posterity, *all* software and versions used when developing the script are provided in the 'development_environment.yml' file located in the 'DOCS' directory for the repository. This file can *not* be used to install the environment on your computer, because many of the dependencies and software builds are system specific. However, this file may help you troubleshoot any dependency errors that may occur in your environment.

#### Installing The Latest Software Versions

1. Create the new environment and specify the software to include in it.
```
conda create -n chipseq-automation -y -c conda-forge -c bioconda \
fastqc bowtie2 samtools macs2 deeptools bedtools sra-tools r-tidyverse r-getopt \
bioconductor-chipseeker bioconductor-rtracklayer bioconductor-genomicranges \
bioconductor-genomicfeatures
```
2. Update the software to the latest compatible versions.
```
conda update -n chipseq-automation -y -c conda-forge -c bioconda --all
```

If you wish to use any of the software in the environment outside of the workflow you can type `conda activate chipseq-automation`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.

## Creating Sample Sheet

In order to keep track of samples, this workflow requires the creation of a sample sheet. An example sheet 'samples.tsv' is provided in the 'DOCS' repository directory. It is important to follow exact formatting of this sheet, as the information within it is used in various stages of the workflow.

| Column | Description |
| ------ | ----------- |
| sample_ID | Short sample identifier (e.g. A001). |
| condition | Experimental condition (e.g. EWSR1_KD). |
| replicate | Sample replicate number (e.g. 1). |
| R1 | Name of R1 fastq file of experimental condition. |
| R2 | Name of R2 fastq file of experimental condition (leave blank if single end). |
| control_ID | Name of input/control sample (e.g. KD_input_1) (leave blank if there is no input/control). |
| R1_control | Name of R1 fastq of input/control (leave blank if there is no input/control). |
| R2_control | Name of R2 fastq of input/control (leave blank if input/control is single end). |

After creating the sample sheet, set the path and file name in the 'settings.conf' file.

## Specifying Run Settings

The last step is to set a few settings in the 'settings.conf' file in the main repository directory. An example settings file is provided in the 'DOCS' directory of the repository.

| Setting | Description |
| ------- | ----------- |
| BASEDIR | The directory for the ChIPseq repository (e.g. /analysis/ChIPseq). |
| OUTDIR | The output directory for the workflow results (e.g. /analysis/results). |
| CORES | The number of CPU cores/threads (e.g. 2). |
| SAMPLE_SHEET | The directory and name of the sample sheet (e.g. /analysis/samples.tsv). |
| DOWNLOAD | Whether the files need to be downloaded from SRA (e.g.'TRUE'). |
| SEQDIR | The directory containing the fastq files (e.g. /analysis/sequences). |
| GENOME_FASTA | The directory and name of the genome assembly fasta (e.g. /analysis/genome/genome.fasta). |
| GENOME_GTF | The directory and name of the genome annotation GTF/GFF (e.g. /analysis/genome/genes.gtf). |
| GENOME_SIZE | Effective genome size (e.g. hg38 is ~3000000000). |
| UPSTREAM | Bases upstream of TSS for defining promoter (e.g. 1000). |
| DOWNSTREAM | Bases downstream of TSS for defining promoter (e.g. 1000). |

## Running the Workflow

After getting the conda environment ready, the sample sheet prepared, and the settings specified, you are now ready to run the workflow. Navigate to the main directory and enter 'bash main.sh'.

###### Notes for IU Folks
If you wish to submit the workflow to a compute node, you can do so by submitting it through the TORQUE resource manager `qsub -l nodes=1:ppn=8,vmem=64gb,walltime=12:00:00 main.sh`. 'ppn' specifies the threads, and 'vmem' is the virtual memory.

# Built With

This workflow would not be possible without the great software listed below.

- [Anaconda](https://www.anaconda.com/) - Software package manager and virtual environment.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Read quality control.
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - Short read aligner.
- [Samtools](http://www.htslib.org/) - SAM/BAM manipulation.
- [MACS](https://github.com/taoliu/MACS) - Peak caller.
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) - Bed file manipulation.
- [deepTools](https://deeptools.readthedocs.io/en/develop/) - Read normalization and heatmap generation.
- [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) - Peak annotation.
- [Tidyverse](https://www.tidyverse.org/) - Data manipulation and visualization in R.
- [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) - Easy genomic file manipulation.
- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) - Robust data structure for genomics data.
- [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) - Working with genomic annotation files.
