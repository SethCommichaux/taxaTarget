#!/bin/sh
#SBATCH --job-name=strict
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=10G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=1


# first need to split strict classifier files into many
# split -l 100000 ../data/strict_classifiers.txt tmp
# then process each in parallel

module load python/3.8.1

python /lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/createDB_scripts/filter_regions.py -s $1

