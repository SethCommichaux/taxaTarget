from Bio import SeqIO
import sys

all_kmers = sys.argv[1]

f=open('tmpermental','w')

for h,i in enumerate(SeqIO.parse(all_kmers,'fasta')):
	if h % 1000000 == 0:
		f.close()
		f = open(str(h)+'.fasta','w')
		f.write(">"+str(i.description)+"\n"+str(i.seq)+"\n")
	else:
		f.write(">"+str(i.description)+"\n"+str(i.seq)+"\n")


