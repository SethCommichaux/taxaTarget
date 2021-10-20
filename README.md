# taxaTarget: A simple supervised learning approach for classifying microeukaryotes in metagenomic data
taxaTarget is a tool for the classification of eukaryotes from metagenomic reads. The input is a fastq file of raw metagenomic sequencing reads. The reads are first mapped to the database with Kaiju to quickly identify reads that are likely sequenced from eukaryotic marker genes. The binned reads are then more sensitively aligned with Diamond to the marker genes. The Diamond output is then provided as input to classifiers, which were trained using the UniRef100 database and a supervised learning approach. Finally, the results are aggregated together and the taxonomic profile of the metagenomic sample is generated.

# Requirements
* Linux-based operating system
* Python3 (version 3.6 or higher; was tested on 3.7)
  * Pandas
  * Numpy
* [Kaiju](https://github.com/bioinformatics-centre/kaiju)
* [Diamond](https://github.com/bbuchfink/diamond)
 
# Installation
```
git clone https://github.com/SethCommichaux/taxaTarget.git
mkdir data
cd data
wget database
```

# Running taxaTarget
Kaiju and Diamond must be in the environment for taxaTarget to run! The input to taxaTarget is an uncompressed fastq file to Kaiju.

```
kaijux -f /path/to/taxaTarget/data/directory/marker_geneDB.fasta.kaiju.fmi -i reads.fastq -z 12 -m 9 | grep "^C" > kaiju
```

Next, a python script is used to extract the reads that mapped to the marker genes with Kaiju

```
python extract_kaiju_reads.py -k kaiju -s reads.fastq -o kaiju.fasta
```

Once the reads are extracted they are input to Diamond to be aligned more sensitively to the marker genes.

```
diamond blastx --sensitive --min-score 55 --db --query kaiju.fasta --threads 12 --outfmt 6 --out kaiju.fasta.diamond
```
