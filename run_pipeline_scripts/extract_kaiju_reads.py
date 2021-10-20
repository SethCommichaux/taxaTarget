'''
This program is used for extracting the sequences of reads from files 
(fastq or fasta) that have been queried against a database with Kaiju.
It also outputs a basic summary file, i.e. read_file_info.txt, with the
number of reads in the query fasta/fastq file and the mean read length.
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file with no unmapped reads")
parser.add_argument("-s", help="read file queried w/ kaiju (must have suffix .fastq or .fasta)")
parser.add_argument("-o", help="name for output fasta file of extracted sequences")
args = parser.parse_args()

# first record all ids of reads that mapped with Kaiju
mapped_reads = {i.strip().split('\t')[1]:0 for i in open(args.k)}

read_count = 0
av_read_len = 0

with open(args.o,'w') as out:
	'''
	Parse fastq file that was queried with Kaiju and output read sequences
	in fasta format if they mapped to database.
	'''
	if args.s.endswith('fastq'): # parse file as fastq if suffix is .fastq
		for h,i in enumerate(open(args.s)):
			if h%4 == 0:
				id = i[1:].split(' ')[0].strip()
			elif h%4 == 1:
				s = i.strip()
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('fasta'): # parse file as fasta if suffix is fasta
		from Bio import SeqIO
		for i in SeqIO.parse(args.s,'fasta'):
			id = str(i.id)
			s = str(i.seq)
			read_count += 1
			av_read_len += len(s)
			if id in mapped_reads:
				out.write(">"+id+"\n"+s+"\n")

# Output number of reads in fastq and mean read length
with open('read_file_info.txt','w') as out2:
	out2.write(str(read_count)+'\n'+str(av_read_len/read_count))

