import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="kaiju output file")
parser.add_argument("-s", help="query fastq/fasta file used by kaiju")
parser.add_argument("-o", help="path name for output fasta file of extracted sequences")
args = parser.parse_args()

mapped_reads = {i.strip().split('\t')[1].split(' ')[0]:0 for i in open(args.k)}

neg = {}

if args.s.endswith('fastq'):
	for h,i in enumerate(open(args.s)):
		if h%4 == 0:
			id = i[1:].split(' ')[0].strip()
		elif h%4 == 1:
			if id in mapped_reads:
				s = i.strip()
				neg[s] = 0
elif args.s.endswith('fasta'):
	for i in SeqIO.parse(args.s,'fasta'):
		if str(i.id) in mapped_reads:
			neg[str(i.seq)] = 0
elif args.s.endswith('pep'):
        for i in SeqIO.parse(args.s,'fasta'):
                if str(i.id) in mapped_reads:
                        neg[str(i.seq)] = 0

with open(args.o,'w') as out:
	for seq in neg.keys():
		out.write(">neg\n"+seq+"\n")
