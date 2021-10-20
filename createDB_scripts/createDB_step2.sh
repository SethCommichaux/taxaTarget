#!/bin/sh
#SBATCH --job-name=diamond
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=20G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12

module load diamond

diamond blastp --min-score 55 --db $2 --query $1 --threads 12 --outfmt 6 --out $1.diamond
