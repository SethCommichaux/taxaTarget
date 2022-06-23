import sys

report = sys.argv[1]

names = {i.strip().split('|')[1].strip().lower():i.strip().split('|')[0].strip() for i in open('/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/data/fullnamelineage.dmp')}
nodes = {i.strip().split('|')[0].strip():i.strip().split('|')[2].strip() for i in open('/lustre/scratch/Seth.Commichaux/Busco_Protist_Pipeline/data/nodes.dmp')}

with open(report+'.modified','w') as out:
	out.write('Lineage\tTaxa\tRank\tRead_count\tBusco_count\n')
	for i in open(report):	
		lineage,rdcount,buscos = i.strip().split('\t')
		taxa = lineage.split(';')[-2]
		taxaID = names[taxa]
		rank = nodes[taxaID]
		out.write('%s\t%s\t%s\t%s\t%s\n' % (lineage,taxa,rank,rdcount,buscos))


