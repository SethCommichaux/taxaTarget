f = '../data/marker_gene_metadata.txt'

r = {}

for i in open(f):
	tmp = i.strip().split('\t')
	MG,species,genus,family = tmp[1],tmp[3],tmp[4],tmp[5]
	if species not in r: r[species] = {MG}
	else: r[species].update([MG])
	if genus not in r: r[genus] = {MG}
	else: r[genus].update([MG])
	if family not in r: r[family] = {MG}
	else: r[family].update([MG])

with open('MGs_per_species_genus_family.txt','w') as out:
	for k,v in r.items():
		out.write(k+'\t'+str(len(v)))
		for z in v:
			out.write('\t'+z)
		out.write('\n')
