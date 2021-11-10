# This script runs taxaTarget as a pipeline using Python3.
# The Python3 must have the Pandas and Numpy packages installed.
# Also, kaiju and diamond must be in the environment

# Import required libraries
#
import os
import sys
import argparse

# Read in commandline arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="Input fastq file", required = True)
parser.add_argument("-dir", help="Path to taxaTarget directory e.g., /usr/home/taxaTarget", required = True)
parser.add_argument("-p", help="Padding to add to thresholds trained with missing data. Range from 0 to 1, with 1 filtering results most aggressively. Default = 0.5",default = '0.5')
parser.add_argument("-t", help="Number of threads to use; default = 12",default = '12')
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

# parse commandline arguments
#
reads_fastq = args.f
padding = args.p
threads = args.t
out = reads_fastq+'.taxaTarget'

if args.dir.endswith('/'):
        taxaTarget = args.dir
else:
        taxaTarget = args.dir+'/'

# Paths to directories, tools, and databases
#
protist_data=taxaTarget+"data/"
run_pipeline=taxaTarget+"run_pipeline_scripts/"
kaijuDB=protist_data+"/marker_geneDB.fasta.kaiju.fmi"
queryDB=protist_data+"/marker_geneDB.fasta"
kaiju=protist_data+"kaiju/bin/kaijux"
diamond=protist_data+"diamond"

# Create and change to output directory
#
os.system('mkdir '+out)

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
print('Beginning analysis. Mapping reads with Kaiju.')
os.system('%s -f %s -i ../%s -z %s -m 9 | grep "^C" > %s/kaiju' % (kaiju,kaijuDB,reads_fastq,threads,out))


# Extract reads that aligned to binning database
#
print('Extracting reads mapped by Kaiju.')
os.system('python %s/extract_kaiju_reads.py -k %s/kaiju -s %s -o %s/kaiju.fasta' % (run_pipeline,out,reads_fastq,out))

# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
print('Aligning reads with Diamond.')
os.system('%s blastx --sensitive --min-score 55 --db %s --query %s/kaiju.fasta --threads %s --outfmt 6 --out %s/kaiju.fasta.diamond' % (diamond,queryDB,out,threads,out))

# Classify reads
#
os.system('python %s/classify_reads_strict.py -d %s/kaiju.fasta.diamond -dir %s/ -p %s' % (run_pipeline,out,protist_data,padding))

