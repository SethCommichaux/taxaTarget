#!/bin/sh
#SBATCH --job-name=busco
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=120G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12


# Load modules and software paths into environment
#
module load python/3.8.1


PWD=`pwd`
createDB=$PWD
data=$PWD"/../data/"
busco=$PWD"/../busco/"
run_pipeline=$PWD"/../run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Retrieve classification thresholds for each protein
#
# python $createDB/strict_threshold_classifiers.py -d $data/diamond_results/ -m $data/marker_gene_metadata.txt


# Collect summary information about the classification power of marker genes
#
python $createDB/classification_power_per_protein_region.py


# Split strict classifier files into many then process each in parallel
# 
mkdir $data/classifiers
cd $data/classifiers
split -l 100000 $data/strict_classifiers.txt tmp


# Filter regions for species level classification power
#
for i in tmp*;
do python $createDB/filter_regions.py -s $i -o $i.classifiers;
done


# concatenate results
#
cat *classifiers > $data/strict_classifiers_filtered.txt
