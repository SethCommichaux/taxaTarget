# taxaTarget: A simple supervised learning approach for classifying microeukaryotes in metagenomic data
taxaTarget is a tool for the classification of eukaryotes from metagenomic reads. The input is a fastq file of raw metagenomic sequencing reads. The reads are first mapped to the database with Kaiju to quickly identify reads that are likely sequenced from eukaryotic marker genes. The binned reads are then more sensitively aligned with Diamond to the marker genes. The Diamond output is then provided as input to classifiers, which were trained using the UniRef100 database and a supervised learning approach. Finally, the results are aggregated together and the taxonomic profile of the metagenomic sample is generated.

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
  -t  The number of parallel threads to use\
