#!/usr/bin/env python

import argparse
import os
import 

parser = argparse.ArgumentParser(description='sample aware read alignment')

parser.add_argument('-o', '--outDir', action='store', nargs=1, default=os.getcwd(), type=str, required=True, help="output directory")
parser.add_argument('-t', '--threads', action='store', nargs=1, default=1, type=int, required=False, help="number of cores")
parser.add_argument('--paired', action='store_true', nargs=0, type=bool, required=False, help="if paired run")

args = parser.parse_args()

with open('samples.txt') as f:
	samples = [x.strip().split('\t') for x in file if not x.startswith('sample_ID')]

if args.paired:
	for sample in samples:
		command = str(
