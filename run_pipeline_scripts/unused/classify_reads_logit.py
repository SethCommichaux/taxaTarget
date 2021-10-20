import argparse
import ast
import math
from collections import Counter

def taxa_mappings():
  taxaID2Lineage = {}
  for i in open(args.f):
      tmp = i.strip().upper().split('|')
      taxaID = tmp[0].strip()
      taxaName = tmp[1].lower().strip()
      lineage = tmp[2].lower().strip()+taxaName
      taxaID2Lineage[taxaID] = lineage
  return taxaID2Lineage

def LCA(lineages):
    lineages = map(lambda x: x.split(';'),lineages)
    new_lineage = [i[0].strip() for i in zip(*lineages) if len(set(i)) == 1 if '' not in i]
    return new_lineage

def logistic_prob(slope,intercept,x):
	if 'na' in [slope,intercept]: return 0
	return 1/(1+(math.exp(-1*(intercept+slope*x))))

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-c", help="marker gene thresholds file")
parser.add_argument("-t", help="probability threshold",default=0.8,type=float)
parser.add_argument("-f", help="NCBI fullnamelineage.dmp file")
args = parser.parse_args()

taxaID2Lineage = taxa_mappings()

print('collected taxonomy dictionary')

classifiers = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	classifiers[MG_region] = {}

print('diamond file first pass')

for i in open(args.c):
	tmp = i.strip().split('\t')
	region = tmp[0]
	if region in classifiers:
		thresholds = ast.literal_eval(tmp[1])
		classifiers[region] = thresholds

print('collected relevant classifiers')

metadata = {}

for i in open(args.m):
        tmp = i.strip().split('\t')
        prot,MG,species,genus,family = tmp[0],tmp[1],tmp[3],tmp[4],tmp[5]
        metadata[prot] = [MG,species,genus,family]

print('collected relevant metadata for marker genes')

final_read_classifications = {}

with open('classified_reads.txt','w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
		if aln_len >= 30:
			region = str(int(20 * math.floor(min(start,end)/20)))
			MG_region = MG+'_'+region
			if classifiers[MG_region] != {}:
				busco,species,genus,family = metadata[MG]
				mean_bitscore = bitscore/aln_len
				sslope,sint = classifiers[MG_region].get('s',['na','na'])
				sprob = logistic_prob(sslope,sint,mean_bitscore)
				gslope,gint = classifiers[MG_region].get('g',['na','na'])
				gprob = logistic_prob(gslope,gint,mean_bitscore)
				fslope,fint = classifiers[MG_region].get('f',['na','na'])
				fprob = logistic_prob(fslope,fint,mean_bitscore)
				probs = [fprob,gprob,sprob]
				if sprob >= args.t:
					lineage = taxaID2Lineage[species]
					out.write(read+'\t'+'species'+'\t'+busco+'\t'+MG_region+'\t'+str(probs)+'\t'+str(mean_bitscore)+'\n')
				elif gprob >= args.t:
					lineage = taxaID2Lineage[genus]
					out.write(read+'\t'+'genus'+'\t'+busco+'\t'+MG_region+'\t'+str(probs)+'\t'+str(mean_bitscore)+'\n')
				elif fprob >= args.t:
					lineage = taxaID2Lineage[family]
					out.write(read+'\t'+'family'+'\t'+busco+'\t'+MG_region+'\t'+str(probs)+'\t'+str(mean_bitscore)+'\n')
				else: continue
				if read not in final_read_classifications:
					final_read_classifications[read] = [[lineage],[busco]]
				else:
					final_read_classifications[read][0].append(lineage)
					final_read_classifications[read][1].append(busco)

sample_taxa = {}

with open('final_read_classifications.txt','w') as out:
        for k,v in final_read_classifications.items():
                lineages,buscos = v[0],v[1]
                lca = LCA(lineages)
                out.write(k+'\t'+';'.join(lca)+'\t'+','.join(list(set(buscos)))+'\n')
                line = ';'.join(lca)
                if line not in sample_taxa:
                    sample_taxa[line] = [[k],buscos]
                else:
                    sample_taxa[line][0].append(k)
                    sample_taxa[line][1] += buscos

print('final read classifications produced')

aggregate_taxa = {}

with open('validation.txt','w') as out:
  for k,v in sample_taxa.items():
    counts,buscos,len_counts = v[0],set(v[1]),len(v[0])
    if len_counts < 2: continue
    if len(buscos) < 2: continue
    if len(buscos) < 10:
        if len_counts > 1.7**len(buscos): continue
    for i in counts:
        out.write(i+'\t'+k+'\n')
    line = ''
    for j in k.split(';'):
        line += j+';'
        if line not in aggregate_taxa:
            aggregate_taxa[line] = [counts,buscos]
        else:
            aggregate_taxa[line][0] += counts
            aggregate_taxa[line][1] = aggregate_taxa[line][1] | buscos

print('aggregated results. writing out to taxonomic report')

with open('Taxonomic_report.txt','w') as out:
	for k,v in sorted(aggregate_taxa.items()):
		counts,buscos = len(v[0]),v[1]
		out.write(k+'\t'+str(counts)+'\t'+str(len(buscos))+'\n')





