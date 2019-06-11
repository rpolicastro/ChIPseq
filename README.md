# ChIPseq
Automation of ChIP-seq Workflow

## About

This repository will use software in a containerized conda virtual environment to process multiple ChIP-seq samples. Before proceeding make sure you have singularity installed, the settings properly set in the 'settings.conf' file, and the sample sheet properly formatted.

# Getting Started

## Cloning Repository

To get started, you must first clone the ChIPseq automation repository. Navigate to a directory you would like to clone the repository to and enter `git clone https://github.com/rpolicastro/ChIPseq.git`.

## Installing Singularity

Singularity containers are self contained 'boxes' that house the software and other files necessary for the workflow. The container itself will automatically be downloaded, but you must have the Singularity software installed to both download and use the container. Please refer to the [documentation](https://www.sylabs.io/docs/) on their website.

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

After getting Singularity installed, the sample sheet prepared, and the settings specified, you are now ready to run the workflow. Navigate to the main directory and enter 'bash main.sh'.

###### Notes for IU Folks
If you wish to submit the workflow to a compute node, you can do so by submitting it through the TORQUE resource manager. Navigate to the directory that contains both your 'settings.conf' and 'main.sh' files. Create a file called **submit_workflow.sh** with the following contents:

```
#!/bin/bash

## Navigate back to directory containing the 'main.sh' and 'settings.conf' file.
cd $PBS_O_WORKDIR

## Load the singularity module on Carbonate/Karst
module load singularity/3.2.0

## Start the workflow
bash main.sh
```

You can now submit the workflow.`qsub -l nodes=1:ppn=8,vmem=64gb,walltime=12:00:00 submit_workflow.sh`. 'ppn' specifies the threads/cores, and 'vmem' is the virtual memory.

# Built With

This workflow would not be possible without the great software listed below.

- [Anaconda](https://www.anaconda.com/) - Software package manager and virtual environment.
- [Singularity](https://www.sylabs.io/docs/) - Containerize sofware and files.
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
