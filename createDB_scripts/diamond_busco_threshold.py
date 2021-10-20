import argparse
from Bio import SeqIO

def taxa_mappings():
  nodes,taxaID2Lineage,names = {},{},{}
  
  for i in open(args.n):
          tmp = i.strip().split('|')
          id,rank = tmp[0].strip(),tmp[2].strip()
          nodes[id] = rank
  
  for i in open(args.t):
          tmp = i.strip().upper().split('|')
          taxaID = tmp[0].strip()
          taxaName = tmp[1].lower().strip()
          lineage = tmp[2].lower().strip()+taxaName
          names[taxaName] = taxaID
          taxaID2Lineage[taxaID] = [x.strip() for x in lineage.split(';')]
  
  return nodes,names,taxaID2Lineage

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="Diamond output aligning protein sequences to busco genes")
parser.add_argument("-e", help="File with 214 marker genes used by eukdetect")
parser.add_argument("-f", help="Fasta file with protein sequences that were aligned")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-l", help="Marker gene length cutoffs file")
parser.add_argument("-b", help="Marker gene bitscore cutoffs file")
parser.add_argument("-n", help="path to nodes taxonomy file")
args = parser.parse_args()

nodes,names,taxID2Lineage = taxa_mappings()

print('Taxonomic mappings found!')

bitscore_cutoffs = {i.strip().split('\t')[0]:float(i.strip().split('\t')[1]) for i in open(args.b)}
length_cutoffs = {i.strip().split('\t')[0]:float(i.strip().split('\t')[3]) for i in open(args.l)}
eukdetect_MGS = {i.strip() for i in open(args.e)}

print('Cutoffs for marker genes determined!')

results = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	query,subject,aln_len,bitscore = tmp[0].split('-')[0],tmp[1].split('_')[0],float(tmp[3]),float(tmp[-1])
	# if subject not in eukdetect_MGS: continue
	if bitscore >= bitscore_cutoffs[subject]:
		if aln_len >= length_cutoffs[subject]:
			results[query] = [subject,'NA','NA','NA','NA','NA','NA']

print('Diamond results processed!')

# uniprot busco_id len species_id genus_id family_id lineage seq

for i in SeqIO.parse(args.f,'fasta'):
	id,d,s,l = str(i.id),str(i.description),str(i.seq),str(len(i.seq))
	if id in results:
		results[id][1] = l
		results[id][6] = s
		taxaID = d.split('TaxID=')[1].split(' ')[0]
		if taxaID in taxID2Lineage:
			lineage = [x.strip() for x in taxID2Lineage[taxaID]]
			results[id][5] = ';'.join(lineage)
			for j in lineage:
				if nodes.get(names.get(j,0)) != 0:
					tid = names.get(j)
					rank = nodes.get(tid)
					if rank == 'species': results[id][2] = tid
					elif rank == 'genus':  
						if j == 'vertebrata':
							if 'rhodophyta' not in lineage:
								continue
						results[id][3] = tid
					elif rank == 'family': results[id][4] = tid

# uniprot busco_id len species_id genus_id family_id lineage seq

with open('marker_geneDB.fasta','w') as out, open('marker_gene_metadata.txt','w') as out2:
	for k,v in results.items():
		if v[2] == 'NA': continue # skip proteins that lack a species label
		out.write(">"+k+"\n"+v[6]+"\n")
		out2.write(k+'\t'+'\t'.join(v)+'\n')
