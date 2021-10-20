# This script runs taxaTarget as a pipeline using Python3.
# The Python3 must have the Pandas and Numpy packages installed.
# Also, kaiju and diamond must be in the environment

# Import required libraries
#
import os
import sys
import argparse

# Set paths to directories and databases
#
createDB="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/createDB_scripts/"
protist_data="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/data/"
run_pipeline="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/run_pipeline_scripts/"
kaijuDB=protist_data+"/marker_geneDB.fasta.kaiju.fmi"
queryDB=protist_data+"/marker_geneDB.fasta"

# Read in commandline arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="Input fastq file", required = True)
parser.add_argument("-p", help="Threshold padding to add for thresholds trained with missing data",default = '0.5')
parser.add_argument("-t", help="Number of threads to use",default = '12')
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

reads_fastq = args.f
padding = args.p
threads = args.t
out = reads_fastq+'.taxaTarget'

# Create and change to output directory
#
os.system('mkdir '+out)
os.chdir(out)

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
os.system('kaijux -f %s -i ../%s -z %s -m 9 | grep "^C" > kaiju' % (kaijuDB,reads_fastq,threads))


# Extract reads that aligned to binning database
#
os.system('python %s/extract_kaiju_reads.py -k kaiju -s ../%s -o kaiju.fasta' % (run_pipeline,reads_fastq))

# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
os.system('diamond blastx --top 0 --sensitive --min-score 55 --db %s --query kaiju.fasta --threads %s --outfmt 6 --out kaiju.fasta.diamond' % (queryDB,threads))

# Classify reads
#
os.system('python %s/classify_reads_strict.py -d kaiju.fasta.diamond -dir %s/ -p %s' % (run_pipeline,protist_data,padding))
