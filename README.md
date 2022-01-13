# taxaTarget: A simple supervised learning approach for classifying microeukaryotes in metagenomic data
taxaTarget is a tool for the classification of eukaryotes from metagenomic reads using a marker gene database and classification thresholds determined through supervised learning. The input is a fastq file of raw metagenomic sequencing reads. The reads are first mapped to the database with Kaiju to quickly identify reads that are likely sequenced from eukaryotic marker genes. The binned reads are then more sensitively aligned with Diamond to the marker genes. The Diamond output is then provided as input to classifiers, which were trained using the UniRef100 database and a supervised learning approach. Finally, the results are aggregated together and the taxonomic profile of the metagenomic sample is generated.

# Publication

# Requirements
* Linux-based operating system
* Python3 (version 3.8 or higher)
  * Pandas
  * Numpy
* Kaiju v1.8.2 (https://github.com/bioinformatics-centre/kaiju)
* Diamond v2.0.13 (https://github.com/bbuchfink/diamond)

To install Kaiju (if there are issues installing check the github documentation),
```
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make
```
To install Diamond,
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.13/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```
In some environments Diamond v2.0.13 doesn't work and the author is working to patch in the next release. If this is the case use v2.0.11. If there are still issues check the Diamond documentation.
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```

# Installation of taxaTarget and database
1) install the taxaTarget scripts.
```
git clone https://github.com/SethCommichaux/taxaTarget.git
cd taxaTarget
```
2) Download and decompress the taxaTarget database (the database is about 3.5 GB so it might take a few minutes).
```
wget https://obj.umiacs.umd.edu/taxatarget/data.zip
unzip data.zip
```
3) Index the taxaTarget database for use by Kaiju and Diamond.
```
cd data/

# kaiju
/path/to/kaiju/bin/kaiju-mkbwt -n 8 -o marker_geneDB.fasta.kaiju marker_geneDB.fasta
/path/to/kaiju/bin/kaiju-mkfmi marker_geneDB.fasta.kaiju
rm marker_geneDB.fasta.kaiju.bwt marker_geneDB.fasta.kaiju.sa

# diamond
/path/to/diamond makedb --in marker_geneDB.fasta --db marker_geneDB.fasta --threads 8
```
4) The full paths to the Kaijux executable (in kaiju/src/kaijux), Diamond executable, and the taxaTarget directory need to be updated in the taxaTarget/run_pipeline_scripts/environment.txt file.

# Running taxaTarget on test data
Once you've installed taxaTarget you can check that everything is working by running it on the test data in run_pipeline_scripts/test_data/. The test data consists of paired end reads sequenced from the genome of Lachancea thermotolerans, which is a budding yeast.
```
cd /path/to/taxaTarget/run_pipeline_scripts/
python run_protist_pipeline_fda.py -r test_data/ERR2886542_1.fastq -r2 test_data/ERR2886542_2.fastq -e environment.txt -t 12 -o test_results
```

# Running taxaTarget
The master script (run_protist_pipeline_fda.py) for running taxaTarget is in the run_pipeline_scripts directory.

Basic usage:
```
python run_protist_pipeline_fda.py -r reads_1.fastq -r2 reads_2.fastq -e environment.txt -o results -t 8
```

Full list of parameters:\
-h, --help &nbsp; Show help message\
-r &nbsp; Read file in .fastq or .fastq.gz format. Use with single-end or forward reads and -r2 for reverse reads if paired end.\
-r2 &nbsp; Reverse reads in .fastq or .fastq.gz format.\
-e &nbsp; Path to environment.txt\
-p &nbsp; The amount of threshold padding to add for thresholds trained with missing data, 0.5 by default\
-o &nbsp; Path for desired output directory. Default uses path to reads with .taxaTarget suffix\
-t &nbsp; The number of parallel threads to use, 12 by default\
--tmp &nbsp; If set, will delete tmp and intermediate files.

# Understanding the output
The taxonomic classification results are output in several files to a directory named after the reads fastq with suffix ".taxaTarget":

Here are the files found in the output directory:

1) classified_reads.txt --> the initial classification results for individual reads that get classified as eukaryotic–this allows users to explore the raw results before filtering.
2) marker_gene_read_counts_per_taxa.txt --> a matrix of read counts per marker gene for all eukaryotic lineages detected in the sample after the results are filtered to build the taxonomic profile.
3) final_read_classifications.txt --> the final classification results for all reads included in the taxonomic profile.
4) Taxonomic_report.txt --> the final aggregate taxonomic profile (i.e. read counts per taxonomic lineage and the number of marker genes with mapped reads) of the sample provided.
5) kaiju --> the kaiju output for reads that mapped to the marker genes
6) kaiju.fasta --> a fasta file of the reads that mapped to the marker genes with kaiju
7) kaiju.fasta.diamond --> diamond results for the reads in kaiju.fasta aligned to the marker genes
8) kaiju.fasta.diamond.filtered --> diamond results filtered for best hits based upon mean bit score
9) read_file_info.txt --> lists number of reads and mean read length for input fastq file
