'''
This program is used for extracting the sequences of reads from files 
(fastq or fasta) that have been queried against a database with Kaiju.
It also outputs a basic summary file, i.e. read_file_info.txt, with the
number of reads in the query fasta/fastq file and the mean read length.
'''

import argparse
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file with no unmapped reads")
parser.add_argument("-s", help="single end or forward read file queried w/ kaiju (.fastq or .fastq.gz)")
parser.add_argument("-s2", help="reverse read file queried w/ kaiju (.fastq or .fastq.gz)")
parser.add_argument("-o", help="path to output directory")
args = parser.parse_args()

# first record all ids of reads that mapped with Kaiju
mapped_reads = {i.strip().split('\t')[1].split('#')[0]:0 for i in open(args.k)}

read_count = 0
av_read_len = 0

with open(args.o+'/kaiju.fasta','w') as out:
	'''
	Parse fastq file that was queried with Kaiju and output read sequences
	in fasta format if they mapped to database.
	'''
	if args.s.endswith('fastq') or args.s.endswith('fq'): # parse file as uncompressed fastq
		for h,i in enumerate(open(args.s)):
			if h%4 == 0:
				id = i[1:].split(' ')[0].strip().split('#')[0]
			elif h%4 == 1:
				s = i.strip()
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
		if args.s2 != None:
			for h,i in enumerate(open(args.s2)):
				if h%4 == 0:
					id = i[1:].split(' ')[0].strip().split('#')[0]
				elif h%4 == 1:
					s = i.strip()
					av_read_len += len(s)
					read_count += 1.0
					if id in mapped_reads:
						out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('fastq.gz') or args.s.endswith('fq.gz'): # parse fastq that is gzipped
		with gzip.open(args.s, "rt") as handle:
			for h,i in enumerate(handle):
				if h%4 == 0:
					id = i[1:].split(' ')[0].strip().split('#')[0]
				elif h%4 == 1:
					s = i.strip()
					av_read_len += len(s)
					read_count += 1.0
					if id in mapped_reads:
						out.write(">"+id+"\n"+s+"\n")
		if args.s2 != None:
			with gzip.open(args.s2, "rt") as handle:
				for h,i in enumerate(handle):
					if h%4 == 0:
						id = i[1:].split(' ')[0].strip().split('#')[0]
					elif h%4 == 1:
						s = i.strip()
						av_read_len += len(s)
						read_count += 1.0
						if id in mapped_reads:
							out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('.fasta') or args.s.endswith('.fa'): # parse file as uncompressed fasta
		for i in SeqIO.parse(args.s,'fasta'):
			id,s = str(i.id).split('#')[0],str(i.seq)
			av_read_len += len(s)
			read_count += 1.0
			if id in mapped_reads:
				out.write(">"+id+"\n"+s+"\n")
		if args.s2 != None:
			for i in SeqIO.parse(args.s2,'fasta'):
				id,s = str(i.id).split('#')[0],str(i.seq)
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('.fasta.gz') or args.s.endswith('.fa.gz'): # parse file as gzip compressed fasta
		with gzip.open(args.s, "rt") as handle:
			for i in SeqIO.parse(handle,'fasta'):
				id,s = str(i.id).split('#')[0],str(i.seq)
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
		if args.s2 != None:
			with gzip.open(args.s2, "rt") as handle:
				for i in SeqIO.parse(handle,'fasta'):
					id,s = str(i.id).split('#')[0],str(i.seq)
					av_read_len += len(s)
					read_count += 1.0
					if id in mapped_reads:
						out.write(">"+id+"\n"+s+"\n")
	else:
		sys.exit('No fasta/fastq file was recognized! Make sure read files have .fastq/.fasta/.fa/.fq suffix. Can also be gzipped with .gz at end.')

# Output number of reads in fastq and mean read length
with open(args.o+'/read_file_info.txt','w') as out2:
	out2.write(str(read_count)+'\n'+str(av_read_len/read_count))

