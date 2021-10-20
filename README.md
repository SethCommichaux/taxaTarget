# taxaTarget: A simple supervised learning approach for classifying microeukaryotes in metagenomic data
taxaTarget is a tool for the classification of eukaryotes from metagenomic reads using a marker gene database and classification thresholds determined through supervised learning. The input is a fastq file of raw metagenomic sequencing reads. The reads are first mapped to the database with Kaiju to quickly identify reads that are likely sequenced from eukaryotic marker genes. The binned reads are then more sensitively aligned with Diamond to the marker genes. The Diamond output is then provided as input to classifiers, which were trained using the UniRef100 database and a supervised learning approach. Finally, the results are aggregated together and the taxonomic profile of the metagenomic sample is generated.

# Publication

# Requirements
* Linux-based operating system
* Python3 (version 3.6 or higher; was tested on 3.8)
  * Pandas
  * Numpy
* [Kaiju](https://github.com/bioinformatics-centre/kaiju)
* [Diamond](https://github.com/bbuchfink/diamond)
 
# Installation
```
git clone https://github.com/SethCommichaux/taxaTarget.git
cd taxaTarget
mkdir data
cd data
wget database
```

# Running taxaTarget
Once installed, you'll find the master script for running taxaTarget in the run_pipeline_scripts directory.

Basic usage:
```
python /path/to/run_pipeline_scripts/run_protist_pipeline_fda.py -f reads.fastq
```

Optionally, there are two parameters that can be set by the user:\
  -p  The amoung of threshold padding to add for thresholds trained with missing data; by default, 0.5\
  -t  The number of parallel threads to use; by default, 12\

Example:
```
python /path/to/run_pipeline_scripts/run_protist_pipeline_fda.py -f reads.fastq -t 8 -p 0.5
```

# Understanding the output
The taxonomic classification results are output in four files:
1) classified_reads.txt provides the initial classification results for individual reads that get classified as eukaryotic–this allows users to explore the raw results before filtering.
2) marker_gene_read_counts_per_taxa.txt is a matrix of read counts per marker gene for all eukaryotic lineages detected in the sample after the results are filtered to build the taxonomic profile.
3) final_read_classifications.txt provides the final classification results for all reads included in the taxonomic profile.
4) Taxonomic_report.txt provides the final aggregate taxonomic profile (i.e. read counts per taxonomic lineage and the number of marker genes with mapped reads) of the sample provided.


