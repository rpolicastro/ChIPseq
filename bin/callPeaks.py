#!/usr/bin/env python

import argparse
import os
import subprocess 
import pandas as pd

### recieve command line arguments

parser = argparse.ArgumentParser(description='sample aware read alignment')

parser.add_argument('--outdir', action='store', nargs=1, default=os.getcwd(), type=str, required=True, help='output directory')
parser.add_argument('--threads', action='store', nargs=1, default=1, type=int, required=False, help='number of cores')
parser.add_argument('--samplesheet', actions='store', nargs=1, type=str, required=True, help='sample sheet file')

args = parser.parse_args()

### functions

def call_peaks(row):
	if row['paired'].lower() == 'paired':
		command = [
			'macs2 callpeak',
			'-t', os.path.join(args.outdir, 'results', 'aligned', row['sample_ID'] + '_' + row['condition'] + '_' + row['replicate'] + '.bam'),
			'-c', os.path.join(args.outdir, 'results', 'aligned', row['control_ID'] + '.bam'),
			'-n', row['sample_ID'] + '_' + row['condition'] + '_' + row['replicate']),
			'-o', os.path.join(args.outdir, 'results', 'peaks'),
			'-f BAMPE -g 3e9'
		]
		subprocess.run(' '.join(command), shell=True, check=True)
	else:
		command = [
			'macs2 callpeak',
			'-t', os.path.join(args.outdir, 'results', 'aligned', row['sample_ID'] + '_' + row['condition'] + '_' + row['replicate'] + '.bam'),
			'-c', os.path.join(args.outdir, 'results', 'aligned', row['control_ID'] + '.bam'),
			'-n', row['sample_ID'] + '_' + row['condition'] + '_' + row['replicate']),
			'-o', os.path.join(args.outdir, 'results', 'peaks'),
			'-f BAM -g 3e9'
		]
		subprocess.run(' '.join(command), shell=True, check=True)

### load sample sheet and call peaks

samples = pd.read_csv(args.samplesheet, sep='\t', header=0, index_col=False)
samples.apply(call_peaks, axis=1)
