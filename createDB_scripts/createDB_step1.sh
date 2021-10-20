#!/bin/sh
#SBATCH --job-name=busco
#SBATCH -t 7-0:00              # time limit: (D-HH:MM)
#SBATCH --mem=120G            # memory per node in MB
#SBATCH --nodes=1              # number of nodes
#SBATCH --cpus-per-task=12


# Load modules and software paths into environment
#
module load python/3.8.1
module load diamond
module load kaiju


PWD=`pwd`
createDB=$PWD
data=$PWD"/../data/"
busco=$PWD"/../busco/"
run_pipeline=$PWD"/../run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Download, decompress NCBI taxonomy files
#
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# tar xvf new_taxdump.tar.gz
# rm names.dmp new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Download BUSCO datasets for eukaryota, fungi and protists
#
# wget https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz
# tar xvf eukaryota_odb10.2020-09-10.tar.gz
# rm *.tar.gz


# Download UniRef100 protein sequence file
#
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz


# Extract eukarota and non-eukaryota protein sequences
#
# python $createDB/splitUniRef100.py -g uniref100.fasta.gz -t fullnamelineage.dmp -n nodes.dmp


# Identify eukaryota marker genes
#
# diamond makedb --in eukaryota.pep --db eukaryota.pep --threads 12
# diamond blastp --threads 12 --query eukaryota.pep --db eukaryota_odb10/ancestral_variants --outfmt 6 --out eukaryota.pep.diamond
# python $createDB/diamond_busco_threshold.py -b $data/eukaryota_odb10/scores_cutoff -l $data/eukaryota_odb10/lengths_cutoff -d $data/eukaryota.pep.diamond -e $createDB/214_eukdetect_MGs.txt -f $data/eukaryota.pep -t $data/fullnamelineage.dmp -n $data/nodes.dmp
# rm eukaryota.pep*


# make marker gene database indexes for diamond and kaiju softwares
#
# diamond makedb --in $data/marker_geneDB.fasta --db $data/marker_geneDB.fasta --threads 12
# mkbwt -o $data/marker_geneDB.fasta.kaiju -n 12 -l 100000 $data/marker_geneDB.fasta
# mkfmi $data/marker_geneDB.fasta.kaiju
# rm $data/marker_geneDB.fasta.kaiju.bwt $data/marker_geneDB.fasta.kaiju.sa


# Identify candidates for non-protist/non-fungi negatives dataset
#
# kaijup -a greedy -e 3 -f $data/marker_geneDB.fasta.kaiju.fmi -i $data/non_eukaryota.pep -z 12 -m 12 -s 50 | grep "^C" > noneuk2euk
# python $createDB/extract_kaiju_reads.py -k noneuk2euk -s non_eukaryota.pep -o negatives.pep


# Extract kmers from positive and negative datasets
#
# python $createDB/kmer_split.py -k 70 -s 10 -i $data/negatives.pep -o $data/negatives70_10.pep
# python $createDB/kmer_split.py -k 70 -s 10 -i $data/marker_geneDB.fasta -o $data/marker_geneDB70_10.pep


# Prepare reads for parallel diamond alignment
#
# mkdir diamond_results
# cat $data/negatives70_10.pep $data/marker_geneDB70_10.pep > $data/diamond_results/all70_10.pep
cd diamond_results
# python $createDB/diamond_preprocess_kmers.py all70_10.pep


# Submit parallel jobs for diamond alignment
#
for i in $data/diamond_results/*.fasta
do
sbatch $createDB/createDB_step2.sh $i $data/marker_geneDB.fasta
done 
