import numpy as np
from sklearn.linear_model import LogisticRegression
import argparse
import ast
import math
from collections import Counter

def LCA(lineages):
    lineages = map(lambda x: x.split(';'),lineages)
    new_lineage = [i[0].strip() for i in zip(*lineages) if len(set(i)) == 1 if '' not in i]
    return new_lineage

def train_logistic_regression(x,labels):
        x = np.array(x)
        X = x[:, np.newaxis]
        clf = LogisticRegression(multi_class='multinomial', solver='newton-cg').fit(X, labels)
        return clf

def taxa_mappings():
  taxaID2Lineage = {}
  for i in open(args.f):
      tmp = i.strip().upper().split('|')
      taxaID = tmp[0].strip()
      taxaName = tmp[1].lower().strip()
      lineage = tmp[2].lower().strip()+taxaName
      taxaID2Lineage[taxaID] = lineage
  return taxaID2Lineage

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-p", help="multilogit data file")
parser.add_argument("-f", help="NCBI fullnamelineage.dmp file")
parser.add_argument("-t", help="confidence score threshold",default=0.2,type=float)
args = parser.parse_args()

taxaID2Lineage = taxa_mappings()

print('collected taxonomy dictionary')

classifiers = {}

for i in open(args.d):
        tmp = i.strip().split('\t')
        read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
        if aln_len < 30: continue
        region = str(int(20 * math.floor(min(start,end)/20)))
        MG_region = MG+'_'+region
        classifiers[MG_region] = {}

print('identified classifiers to be used for analysis')

metadata = {}

for i in open(args.m):
	tmp = i.strip().split('\t')
	prot,MG,species,genus,family = tmp[0],tmp[1],tmp[3],tmp[4],tmp[5]
	metadata[prot] = [MG,species,genus,family]

print('collected metadata for classifiers')

for i in open(args.p):
	tmp = i.strip().split('\t')
	region,x,labels,dt = tmp[0],[float(j) for j in tmp[1].split(',')],tmp[2].split(','),float(tmp[3])
	if region in classifiers:
		classifiers[region] = [x,labels]

print('collected training data...begin processing reads')

trainers = {}
final_read_classifications = {}

with open('classified_reads.txt','w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
		if aln_len < 30: continue
		region = str(int(20 * math.floor(min(start,end)/20)))
		MG_region = MG+'_'+region
		if classifiers.get(MG_region,0) in [0,{}]: continue
		mean_bitscore = bitscore/aln_len
		if MG_region not in trainers:
			trainers[MG_region] = train_logistic_regression(classifiers[MG_region][0],classifiers[MG_region][1])
		probs = trainers[MG_region].predict_proba([[mean_bitscore]])[0]
		classes = trainers[MG_region].classes_
		predicted_class = trainers[MG_region].predict([[mean_bitscore]])[0]
		confidence = trainers[MG_region].decision_function([[mean_bitscore]])[0]
		if type(confidence) != np.float64: confidence = max(confidence)
		if predicted_class != 'n':
			if confidence >= float(args.t):
				busco,species,genus,family = metadata[MG]
				if 'f' == predicted_class: lineage = taxaID2Lineage[family]
				elif 'g' == predicted_class: lineage = taxaID2Lineage[genus]
				elif 's' == predicted_class: lineage = taxaID2Lineage[species]
				out.write(read+'\t'+MG_region+'\t'+busco+'\t'+predicted_class+'\t'+lineage+'\t'+str(probs)+'\t'+str(classes)+'\t'+str(mean_bitscore)+'\n')
				if read not in final_read_classifications:
					final_read_classifications[read] = [lineage]
				else:
					final_read_classifications[read].append(lineage)



'''
		results = dict(zip(classes, probs))
		best_class,best_prob = 0,0
		for r in ['s','g','f','n']:
			if results.get(r,0) > best_prob:
				best_prob = results[r]
				best_class = r
		if best_class != 'n':
			second_best = sorted(probs)[-2]
			if (best_prob - second_best) >= 0.2/len(probs):
				busco,species,genus,family = metadata[MG]
				if 'f' == best_class:
					out.write(read+'\t'+MG_region+'\t'+busco+'\t'+best_class+'\t'+taxaID2Lineage[family]+'\t'+str(probs)+'\t'+str(classes)+'\t'+str(mean_bitscore)+'\n')
					lineage = taxaID2Lineage[family]
				elif 'g' == best_class:
					out.write(read+'\t'+MG_region+'\t'+busco+'\t'+best_class+'\t'+taxaID2Lineage[genus]+'\t'+str(probs)+'\t'+str(classes)+'\t'+str(mean_bitscore)+'\n')
					lineage = taxaID2Lineage[genus]
				elif 's' == best_class:
					out.write(read+'\t'+MG_region+'\t'+busco+'\t'+best_class+'\t'+taxaID2Lineage[species]+'\t'+str(probs)+'\t'+str(classes)+'\t'+str(mean_bitscore)+'\n')
					lineage = taxaID2Lineage[species]
'''

sample_taxa = Counter([])

with open('final_read_classifications.txt','w') as out:
	for k,v in final_read_classifications.items():
		lca = LCA(v)
		out.write(k+'\t'+';'.join(lca)+'\n')
		line = ''
		for j in lca:
			line += j+';'
			sample_taxa.update([line])

with open('Taxonomic_report.txt','w') as out:
	for k,v in sorted(sample_taxa.items()):
		out.write(k+'\t'+str(v)+'\n')

