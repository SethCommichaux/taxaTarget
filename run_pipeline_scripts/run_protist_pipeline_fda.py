# This script runs taxaTarget as a pipeline using Python3.
# The Python3 must have the Pandas and Numpy packages installed.
# Also, kaiju and diamond must be in the environment

# check for required libraries in python environment
#
import os
import sys
import argparse
import math
import gzip
import random
from collections import Counter

try:
	import pandas as pd
except ImportError:
	sys.exit('Pandas library not detected in Python environment!')
try:
	import numpy as np
except ImportError:
	sys.exit('Numpy library not detected in Python environment!')

# Read in commandline arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("-r", help="Read file in .fastq or .fastq.gz format. Use with single-end or forward reads and -r2 for reverse reads if paired end", required = True)
parser.add_argument("-r2", help="Reverse reads in .fastq or .fastq.gz format")
parser.add_argument("-e", help="Path to environment.txt", required = True)
parser.add_argument("-p", help="Padding to add to thresholds trained with missing data. Range from 0 to 1, with 1 filtering results most aggressively. Default = 0.5",default = '0.5')
parser.add_argument("-o", help="Path for desired output directory. Default uses path to reads with .taxaTarget suffix")
parser.add_argument("-t", help="Number of threads to use; default = 12",default = '12')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

# parse commandline arguments
#
reads_fastq = args.r
reverse_reads = args.r2

padding = args.p
threads = args.t

environment = {}

for i in open(args.e):
	if i[0] != '#':
		k,v = i.strip().replace("'","").split('=')
		environment[k] = v 

if args.o == None:
	out = reads_fastq+'.taxaTarget'
else:
	out = args.o

# Paths to directories, tools, and databases
#
protist_data=environment['taxatarget']+"data/"
run_pipeline=environment['taxatarget']+"run_pipeline_scripts/"
kaijuDB=protist_data+"/marker_geneDB.fasta.kaiju.fmi"
queryDB=protist_data+"/marker_geneDB.fasta"
kaiju=environment['kaiju']
diamond=environment['diamond']

# Create and change to output directory
#
os.system('mkdir '+out)

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
print('Beginning analysis. Mapping reads with Kaiju.')
if reverse_reads == None:
	os.system('%s -f %s -i %s -z %s -m 9 | grep "^C" > %s/kaiju' % (kaiju,kaijuDB,reads_fastq,threads,out))
else:
	os.system('%s -f %s -i %s -j %s -z %s -m 9 | grep "^C" > %s/kaiju' % (kaiju,kaijuDB,reads_fastq,reverse_reads,threads,out))

# Extract reads that aligned to binning database
#
print('Extracting reads mapped by Kaiju.')
if reverse_reads == None:
	os.system('python %s/extract_kaiju_reads.py -k %s/kaiju -s %s -o %s' % (run_pipeline,out,reads_fastq,out))
else:
	os.system('python %s/extract_kaiju_reads.py -k %s/kaiju -s %s -s2 %s -o %s' % (run_pipeline,out,reads_fastq,reverse_reads,out))

# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
print('Aligning reads with Diamond.')
os.system('%s blastx --sensitive --min-score 55 --db %s --query %s/kaiju.fasta --threads %s --outfmt 6 --out %s/kaiju.fasta.diamond' % (diamond,queryDB,out,threads,out))

# Classify reads
#
os.system('python %s/classify_reads_strict.py -d %s/kaiju.fasta.diamond -dir %s/ -p %s -o %s' % (run_pipeline,out,protist_data,padding,out))
