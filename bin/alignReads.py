#!/usr/bin/env python

import argparse
import os
import subprocess 
import pandas as pd

### recieve command line arguments

parser = argparse.ArgumentParser(description='sample aware read alignment')

parser.add_argument('--outDir', action='store', nargs=1, default=os.getcwd(), type=str, required=True, help='output directory')
parser.add_argument('--threads', action='store', nargs=1, default=1, type=int, required=False, help='number of cores')
parser.add_argument('--sampleSheet', actions='store', nargs=1, type=str, required=True, help='sample sheet file')

args = parser.parse_args()

### functions

def align_experimental(row):
	if row['paired'].lower() == 'paired':
		command = [
			'bowtie2',
			'-x', os.path.join(args.outdir, 'genome', 'hg38'),
			'-1', row['R1'], '-2', row['R2'],
			'-S', os.path.join(args.outdir, 'results', 'aligned', row['sampleID'] + '_' + row['condition'] + '_' + row['replicate'] + '.sam'),
			'-q --phred33 --no-mixed --no-discordant --threads', args.threads
		]
		subprocess.run(' '.join(command), shell=True, check=True)
	else:
		command = [
                        'bowtie2',
                        '-x', os.path.join(args.outdir, 'genome', 'hg38'),
                        '-U', row['R1'],
                        '-S', os.path.join(args.outdir, 'results', 'aligned', row['sampleID'] + '_' + row['condition'] + '_' + row['replicate'] + '.sam'),
                        '-q --phred33 --threads', args.threads
                ]
                subprocess.run(' '.join(command), shell=True, check=True)

def align_controls(row):
        if row['paired'].lower() == 'paired':
                command = [
                        'bowtie2',
                        '-x', os.path.join(args.outdir, 'genome', 'hg38'),
                        '-1', row['R1_control'], '-2', row['R2_control'],
                        '-S', os.path.join(args.outdir, 'results', 'aligned', row['control_ID'] + '.sam'),
                        '-q --phred33 --no-mixed --no-discordant --threads', args.threads
                ]
                subprocess.run(' '.join(command), shell=True, check=True)
        else:
                command = [
                        'bowtie2',
                        '-x', os.path.join(args.outdir, 'genome', 'hg38'),
                        '-U', row['R1_control'],
                        '-S', os.path.join(args.outdir, 'results', 'aligned', row['control_ID'] + '.sam'),
                        '-q --phred33 --threads', args.threads
                ]
                subprocess.run(' '.join(command), shell=True, check=True)

### load sample sheet as pandas object

samples = pd.read_csv(args.sampleSheet, sep='\t', header=0, index_col=False)
samples.apply(align_experimental, axis=1)
unique_controls = samples.drop_duplicates('control_ID', keep='first')
unique_controls.apply(align_controls, axis=1)
