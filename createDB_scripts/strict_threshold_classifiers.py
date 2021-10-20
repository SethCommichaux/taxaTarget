import argparse
import os
import math

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="Directory with diamond alignment files --> reads to marker genes mappings")
parser.add_argument("-m", help="Marker gene metadata file")
args = parser.parse_args()

taxa_mappings = {}
results = {}
c = 0

for i in open(args.m):
	tmp = i.strip().split('\t')
	prot,species,genus,family = tmp[0],tmp[3],tmp[4],tmp[5]
	taxa_mappings[prot] = [species,genus,family]

# structure of results
# results[marker_gene_region] = [negative_max,species_min,species_max,genus_min,genus_max,family_min,family_max]

for f in os.listdir(args.d):
	if f.endswith('.diamond'):
		c += 1
		print(c)
		for i in open(args.d+'/'+f):
			tmp = i.strip().split('\t')
			read,MG,aln_len,start,end,bitscore = '_'.join(tmp[0].split('_')[:-1]),tmp[1],float(tmp[3]),float(tmp[8]),float(tmp[9]),float(tmp[-1])
			if aln_len < 30: continue
			region = str(int(20 * math.floor(min(start,end)/20))) # floor round to nearest 20
			MG_region = MG+'_'+region
			if MG_region not in results:
				results[MG_region] = [0,100,0,100,0,100,0]
			mean_bitscore = bitscore/aln_len
			if read == 'neg':
				results[MG_region][0] = max(mean_bitscore,results[MG_region][0])
			else:
				read_lineage = taxa_mappings[read]
				MG_lineage = taxa_mappings[MG]
				for j in range(3):
					if (read_lineage[j] == MG_lineage[j]) and (read_lineage[j] != 'NA'):
						relation = ['species','genus','family'][j]
						if relation == 'species':
							results[MG_region][1] = min(mean_bitscore,results[MG_region][1])
							results[MG_region][2] = max(mean_bitscore,results[MG_region][2])
						elif relation == 'genus':
							results[MG_region][3] = min(mean_bitscore,results[MG_region][3])
							results[MG_region][4] = max(mean_bitscore,results[MG_region][4])
						elif relation == 'family':
							results[MG_region][5] = min(mean_bitscore,results[MG_region][5])
							results[MG_region][6] = max(mean_bitscore,results[MG_region][6])
						break
				else:
					results[MG_region][0] = max(mean_bitscore,results[MG_region][0])


with open('strict_classifiers.txt','w') as out:
	for k,v in results.items():
		v = map(str,v)
		out.write(k+'\t'+'\t'.join(v)+'\n')
