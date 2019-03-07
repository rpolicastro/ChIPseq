# ChIPseq
Automation of ChIP-seq Workflow

## About

This repository will use software in a conda virtual environment to process multiple ChIP-seq samples. Before preceding make sure you have your conda environment installed, and the settings properly set in the 'settings.conf' file.

## Preparing Conda Environment

This workflow takes advantage of the [conda](https://conda.io/en/latest/) package manager and virtual environment. The conda package manager installs both the main software and all dependencies into a 'virtual environment' to ensure compatabilty. Furthermore, the provided 'environment.yml' file will reproduce the software environment used when developing the workflow. This ensures prolonged compatabilty and reproducibility.

Before creating the environment, you must first install miniconda on your computer.
1. [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda), and make sure that conda is in your PATH.
2. Update conda to the latest version `conda update conda`.

You are now read to create the virtual sofware environment, and download all software and dependencies.

1. Create the environment and specify the software to include in it `conda create -n chipseq-automation -y -c conda-forge -c bioconda fastqc bowtie2 samtools macs2 pandas`.
2. Update the software to the latest compatible versions `conda update -n chipseq-automation -y -c conda-forge -c bioconda --all`.

If you wish to use any of the software in the environment outside of the workflow you can type `conda activate chipseq-automation`. You can deactivate the environment by closing your terminal or entering `conda deactivate`.
