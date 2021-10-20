from collections import Counter

f = 'strict_classifiers.txt'
uniprot2MG = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open('../data/marker_gene_metadata.txt')}

r = {}
total = {'s':0,'g':0,'f':0,'n':0}
perMG = {}
p = ['n','s','s','g','g','f','f']

for h,i in enumerate(open(f)):
	if h%10000 == 0: print(h)
	tmp = i.strip().split('\t')
	tmp2 = []
	prot = '_'.join(tmp[0].split('_')[:2])
	if prot not in r:
		r[prot] = []
		perMG[uniprot2MG[prot]] = {'s':0,'g':0,'f':0,'n':0}
	scores = []
	for j in tmp[1:]:
		if j == '100': scores.append(0)
		else: scores.append(float(j))
	m = max(scores)
	for j,k in enumerate(scores):
		if k == m:
			tmp2.append(p[j])
	for j in ['n','f','g','s']:
		if j in tmp2:
			r[prot].append(j)
			total[j] += 1
			perMG[uniprot2MG[prot]][j] += 1
			break


# for k,v in perMG.items():
#	print k,'\t',v['n'],'\t',v['s'],'\t',v['g'],'\t',v['f']

for k,v in total.items():
	print(k,v)

with open('per_protein_region_classification.txt','w') as out:
	out.write('UniProtID\tBUSCO_ID\tNegative\tFamily\tGenus\tSpecies\n')
	for k,v in r.items():
		c = dict(Counter(v))
		mg = uniprot2MG[k]
		out.write(k+'\t'+mg+'\t'+str(c.get('n',0))+'\t'+str(c.get('f',0))+'\t'+str(c.get('g',0))+'\t'+str(c.get('s',0))+'\n')

		

