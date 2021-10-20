import argparse
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to UniRef format protein fasta file that is gzipped")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-n", help="path to nodes taxonomy file")
args = parser.parse_args()


# parse NCBI taxonomy fullnamelineage.dmp file for taxaIDs, taxaNames and taxonomic lineages
def build_taxa_dict():
	taxaID2_Name_Lineage,names,nodes = {},{},{}
	
	for i in open(args.t):
		tmp = [j.strip() for j in i.upper().split('|')]
		taxaID = tmp[0]
		taxaName = tmp[1]
		lineage = tmp[2]
		taxaID2_Name_Lineage[taxaID] = lineage+taxaName
		names[taxaName] = taxaID
	print('processed full lineages!')
	for i in open(args.n):
		tmp = [j.strip() for j in i.split('|')]
		id,rank = tmp[0],tmp[2]
		if rank == 'species':
			nodes[id] = rank
	print('processed nodes!')
	return taxaID2_Name_Lineage,names,nodes


taxaID2_Name_Lineage,names,nodes = build_taxa_dict()

with open('eukaryota.pep','w') as out1, open('non_eukaryota.pep','w') as out2:
	g = gzip.open(args.g,'rt')
	for i in SeqIO.parse(g,'fasta'):
		d = str(i.description)
		s = str(i.seq)
		taxaID = d.split('TaxID=')[1].split(' ')[0]
		if taxaID2_Name_Lineage.get(taxaID,0) != 0:
			lineage = [j.strip() for j in taxaID2_Name_Lineage[taxaID].split(';')]
			if 'EUKARYOTA' in lineage:
				if 'ENVIRONMENTAL SAMPLES' in lineage: continue
				if 'UNCLASSIFIED EUKARYOTES' in lineage: continue
				for group in lineage:
					if nodes.get(names.get(group,0),0) == 'species':
						out1.write(">"+d+"\n"+s+"\n")
						break
			else:
				out2.write(">"+d+"\n"+s+"\n")


