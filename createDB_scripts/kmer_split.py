from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-k", help="Desired kmer size for output queries.pep file")
parser.add_argument("-i", help="Input pep file")
parser.add_argument("-s", help="Desired step size in amino acids")
parser.add_argument("-o", help="Desired name for output file")
args = parser.parse_args()

kmer_len = int(args.k)
step = int(args.s)

kmers = {}

if args.i == 'negatives.pep':
	neg = {}
	for i in SeqIO.parse(args.i,'fasta'):
		id = str(i.id)
		seq = str(i.seq)
		seq_len = len(seq)
		for j in range(0,seq_len+1-step,step):
			kmer = seq[j:j+kmer_len]
			if len(kmer) < kmer_len: continue
			if abs(j-seq_len) < kmer_len + step: # e.g. abs(323 - 400) < 70 + 10 == 77 < 80
				neg[seq[j:]] = 0
			else:
				neg[kmer] = 0
	with open(args.o,'w') as out:
		for k in neg.keys():
			out.write(">neg\n"+k+"\n")
else:
	with open(args.o,'w') as out:
		for i in SeqIO.parse(args.i,'fasta'):
			id = str(i.id)
			seq = str(i.seq)
			seq_len = len(seq)
			for j in range(0,seq_len+1-step,step):
				kmer = seq[j:j+kmer_len]
				if len(kmer) < kmer_len: continue
				if abs(j-seq_len) < kmer_len + step: # e.g. abs(323 - 400) < 70 + 10 == 77 < 80
					out.write(">"+id+"_"+str(j)+"\n"+seq[j:]+"\n")
				else:
					out.write(">"+id+"_"+str(j)+"\n"+kmer+"\n")

