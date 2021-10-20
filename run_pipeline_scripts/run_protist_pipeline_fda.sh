#!/bin/sh
#SBATCH --job-name=taxaTarget
#SBATCH -t 1-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=20G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12

# Load modules and software paths into environment
#
module load diamond
module load kaiju
module load python/3.8.1


# Set paths to directories and databases
#
createDB="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/createDB_scripts/"
protist_data="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/data/"
run_pipeline="/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/run_pipeline_scripts/"
kaijuDB=$protist_data"/marker_geneDB.fasta.kaiju.fmi"
queryDB=$protist_data"/marker_geneDB.fasta"


################################################################
################################################################
################################################################


echo taxaTarget $1

# Fastq file(s) to be analyzed
#
reads_fastq=$1

# Output file
#
out=$1.taxaTarget
mkdir $out
cd $out

# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
kaijux -f $kaijuDB -i ../$reads_fastq -z 12 -m 9 | grep "^C" > kaiju


# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k kaiju -s ../$reads_fastq -o kaiju.fasta


# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
diamond blastx --top 0 --sensitive --min-score 55 --db $queryDB --query kaiju.fasta --threads 12 --outfmt 6 --out kaiju.fasta.diamond


# Classify reads
#
python $run_pipeline/classify_reads_strict.py -d kaiju.fasta.diamond -dir $protist_data/ -p 0.5
