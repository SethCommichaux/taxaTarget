'''
This program is used for extracting the sequences of reads from files 
(fastq or fasta) that have been queried against a database with Kaiju.
It also outputs a basic summary file, i.e. read_file_info.txt, with the
number of reads in the query fasta/fastq file and the mean read length.
'''

import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file with no unmapped reads")
parser.add_argument("-s", help="single end or forward read file queried w/ kaiju (.fastq or .fastq.gz)")
parser.add_argument("-s2", help="reverse read file queried w/ kaiju (.fastq or .fastq.gz)")
parser.add_argument("-o", help="path to output directory")
args = parser.parse_args()

# first record all ids of reads that mapped with Kaiju
mapped_reads = {i.strip().split('\t')[1]:0 for i in open(args.k)}

read_count = 0
av_read_len = 0

with open(args.o+'/kaiju.fasta','w') as out:
	'''
	Parse fastq file that was queried with Kaiju and output read sequences
	in fasta format if they mapped to database.
	'''
	if args.s.endswith('fastq') or args.s.endswith('fq'): # parse file as fastq if suffix is .fastq
		for h,i in enumerate(open(args.s)):
			if h%4 == 0:
				id = i[1:].split(' ')[0].strip()
			elif h%4 == 1:
				s = i.strip()
				av_read_len += len(s)
				read_count += 1.0
				if id in mapped_reads:
					out.write(">"+id+"\n"+s+"\n")
		if args.s2 != None:
			for h,i in enumerate(open(args.s)):
				if h%4 == 0:
					id = i[1:].split(' ')[0].strip()
				elif h%4 == 1:
					s = i.strip()
					av_read_len += len(s)
					read_count += 1.0
					if id in mapped_reads:
						out.write(">"+id+"\n"+s+"\n")
	elif args.s.endswith('fastq.gz') or args.s.endswith('fq.gz'): # parse fastq that is gzipped
                for h,i in enumerate(gzip.open(args.s)):
                        if h%4 == 0:
                                id = i[1:].split(' ')[0].strip()
                        elif h%4 == 1:
                                s = i.strip()
                                av_read_len += len(s)
                                read_count += 1.0
                                if id in mapped_reads:
                                        out.write(">"+id+"\n"+s+"\n")
                if args.s2 != None:
                        for h,i in enumerate(gzip.open(args.s)):
                                if h%4 == 0:
                                        id = i[1:].split(' ')[0].strip()
                                elif h%4 == 1:
                                        s = i.strip()
                                        av_read_len += len(s)
                                        read_count += 1.0
                                        if id in mapped_reads:
                                                out.write(">"+id+"\n"+s+"\n")
	else:
		sys.exit('No fastq file was recognized!')

# Output number of reads in fastq and mean read length
with open(args.o+'/read_file_info.txt','w') as out2:
	out2.write(str(read_count)+'\n'+str(av_read_len/read_count))

