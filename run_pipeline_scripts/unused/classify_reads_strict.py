import argparse
import sys
import math
import pandas as pd
import numpy as np
from collections import Counter

def read_file_info():
	# function extracts number of reads and mean read length from read_file_info.txt
	for h,i in enumerate(open('read_file_info.txt')):
		if h == 0: number_reads = float(i.strip())
		elif h == 1: read_len = float(i.strip())
	return number_reads,read_len

def taxa_mappings():
  # function for parsing NCBI fullnamelineage.dmp and returning a dict mapping taxa ID to taxonomic lineage
  taxaID2Lineage = {}
  for i in open(args.dir+'/fullnamelineage.dmp'):
      tmp = i.strip().upper().split('|')
      taxaID = tmp[0].strip()
      taxaName = tmp[1].lower().strip()
      lineage = tmp[2].lower().strip().replace('; ',';')+taxaName
      if 'eukaryota' in lineage:
          taxaID2Lineage[taxaID] = lineage
  return taxaID2Lineage

def used_classifiers():
	# identifies all 20 amino acid windows of all database genes that will be used for taxonomic classification
	# it reduces memory usage to only load classifier data that is relevant to the sample rather than all classifier data
	c = {}
	for i in open(args.d):
		tmp = i.strip().split('\t')
		MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
		region = str(int(20 * math.floor(min(start,end)/20)))
		MG_region = MG+'_'+region
		c[MG_region] = {}
	return c

def collect_thresholds(classifiers):
	# for every 20 amino acid region to be used for classification, this function collects the thresholds
	# for negatives, family, genus and species levels
	for i in open(args.dir+'/strict_classifiers_filtered.txt'):
		tmp = i.strip().split('\t')
		region = tmp[0]
		if region in classifiers:
			thresholds = []
			for x in tmp[1:]:
				if x == 'na': thresholds.append('na')
				else: thresholds.append(float(x))
			classifiers[region] = thresholds
	return classifiers

def collect_classifier_metadata(taxaID2Lineage):
	# function collects metadata associated with each classifier
	# example metadata includes: UniProt ID, BUSCO marker gene ID, protein length, taxonomic IDs for
	# family, genus and species, as well as the full taxonomic lineage
	metadata,mg_per_species,prot2len = {},{},{}
	
	for i in open(args.dir+'/marker_gene_metadata.txt'):
		tmp = i.strip().split('\t')
		prot,MG,prot_len,species,genus,family,lineage = tmp[0],tmp[1],float(tmp[2]),tmp[3],tmp[4],tmp[5],tmp[6]
		metadata[prot] = [MG,species,genus,family]
		species_lineage = taxaID2Lineage[species]
		prot2len[prot] = prot_len
		if species_lineage not in mg_per_species: mg_per_species[species_lineage] = set([MG])
		else: mg_per_species[species_lineage].update([MG])
	return metadata,mg_per_species,prot2len

def pad_thresholds(nmax,smin,gmin,fmin,padding):
	if nmax == 0:
		# if region lacked negative training data, pad family, genus and species thresholds
		smin += padding
		if type(gmin) != str: gmin += padding
		if type(fmin) != str: fmin += padding
	elif fmin == 'na':
		# if region lacked family training data, pad genus and species thesholds 
		smin += padding
		if type(gmin) != str: gmin += padding
	elif gmin == 'na':
		# if region lacked genus training data, pad species threshold
		smin += padding
	return(nmax,smin,gmin,fmin)

def classify_read_first_pass(smin,gmin,fmin,mean_bitscore):
	lineage = ''
	rank = ''
	if mean_bitscore >= smin:
		lineage = taxaID2Lineage[species]
		rank = 'species'
		if genus not in genus_species:
			genus_species[genus] = set([lineage])
		else:
			genus_species[genus].update([lineage])
	elif gmin != 'na':
		if mean_bitscore >= gmin:
			lineage = taxaID2Lineage[genus]
			rank = 'genus'
	elif fmin != 'na':
		if mean_bitscore >= fmin:
			lineage = taxaID2Lineage[family]
			rank = 'family'
	return(lineage,rank)

def LCA(lineages):
    # function for finding the lowest common ancestor of a list of lineages
    lineages = [x.split(';') for x in lineages]
    new_lineage = [i[0].strip() for i in zip(*lineages) if len(set(i)) == 1 if '' not in i]
    return new_lineage

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file", required = True)
parser.add_argument("-dir", help="path to data directory of taxaTarget", required = True)
parser.add_argument("-p", help="proportion of range(0,smax) to add to thresholds for species, genus and family",default=0.02,type=float)
args = parser.parse_args()

# extract number of reads in sample and mean read length
number_reads,read_len = read_file_info()
print('Extracted number of reads and mean read length for sample')

taxaID2Lineage = taxa_mappings()
print('collected taxonomy dictionary')

# maps number of marker genes per family, genus, and species per taxaID
MGs_per_taxa = {taxaID2Lineage[i.strip().split('\t')[0]]:i.strip().split('\t')[2:] for i in open(args.dir+'/MGs_per_species_genus_family.txt')}

classifiers = used_classifiers()
print('diamond file first pass')

classifiers = collect_thresholds(classifiers)
print('collected relevant classifiers')

metadata,MG_per_species,prot2len = collect_classifier_metadata(taxaID2Lineage)
print('collected relevant metadata for marker genes')

read_classifications = {}
multimapped = Counter([i.strip().split('\t')[0] for i in open(args.d)])
genus_species = {}

with open('classified_reads.txt','w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,MG,pident,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[2]),float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
		mean_bitscore = bitscore/aln_len
		if aln_len >= 30:
			region = str(int(20 * math.floor(min(start,end)/20)))
			MG_region = MG+'_'+region
			busco,species,genus,family = metadata[MG]
			if (multimapped[read] > 1) and (pident == 100):
				# assign lowest common ancestor to reads that map ambiguously (i.e. 100% identity) to multiple proteins
				if read not in read_classifications:
					read_classifications[read] = [[taxaID2Lineage[species]],set([busco]),mean_bitscore,[prot2len[MG]]]
				else:
					read_classifications[read][0].append(taxaID2Lineage[species])
					read_classifications[read][1].update([busco])
					read_classifications[read][3].append(prot2len[MG])
			elif classifiers[MG_region] != {}: # check reads map to marker gene regions with classifiers
				nmax,smin,smax,gmin,gmax,fmin,fmax = classifiers[MG_region]
				padding = args.p*smax
				nmax,smin,gmin,fmin = pad_thresholds(nmax,smin,gmin,fmin,padding)
				lineage,rank = classify_read_first_pass(smin,gmin,fmin,mean_bitscore)
				if lineage != '':
					out.write(read+'\t'+rank+'\t'+busco+'\t'+MG_region+'\t'+str(mean_bitscore)+'\t'+lineage+'\n')
					if read not in read_classifications:
						read_classifications[read] = [[lineage],set([busco]),mean_bitscore,[prot2len[MG]]]
					else:
						read_classifications[read][0].append(lineage)
						read_classifications[read][1].update([busco])
						read_classifications[read][3].append(prot2len[MG])

mg_names = [i.strip() for i in open(args.dir+'/255_MGs.txt')] + ['Total_reads']
sample_taxa = {'names':mg_names}
final_read_classifications = {}

for read,v in read_classifications.items():
    lineages,buscos,mean_bitscore,average_prot_len = v[0],v[1],v[2],sum(v[3])/float(len(v[3]))
    lca = LCA(lineages)
    line = ';'.join(lca)
    if line not in sample_taxa:
        sample_taxa[line] = np.zeros(len(mg_names))
    if line not in final_read_classifications:
        final_read_classifications[line] = [read]
    else:
        final_read_classifications[line].append(read)
    for i in buscos:
        sample_taxa[line][mg_names.index(i)] += 1
        sample_taxa[line][-1] += 1 

del_keys = []

for k,v in sample_taxa.items():
    if k == 'names': continue
    lineage_MGS = MGs_per_taxa.get(k,0)
    if lineage_MGS == 0:
        del_keys.append(k)
        del final_read_classifications[k]
        continue
    tmp = []
    non_zero = 0
    for i in lineage_MGS:
        rd_cnt = v[mg_names.index(i)]
        tmp.append(rd_cnt)
        if rd_cnt > 0: non_zero += 1
    total_rd_cnt = sum(tmp)
    expected_number_MGs = 0.96270 * total_rd_cnt + -0.00128 * total_rd_cnt**2 + 1.89510
    if total_rd_cnt < 3:
        print(k)
        del_keys.append(k)
        del final_read_classifications[k]
    elif non_zero < 3:
        del_keys.append(k)
        del final_read_classifications[k]
    elif non_zero < 20:
        if non_zero/float(len(tmp)) < 0.25: # if less than 25% of taxa's marker gene set have mapped reads
            if non_zero < expected_number_MGs*0.25: # if less than 25% of expected number of marker genes have mapped reads
                del_keys.append(k)
                del final_read_classifications[k]

for i in del_keys:
	del sample_taxa[i]

df = pd.DataFrame(sample_taxa)
t = df.transpose()

t.to_csv('marker_gene_read_counts_per_taxa.txt',sep='\t')

print('marker gene counts collected')

for k,v in genus_species.items():
	s = {}
	for line in v:
		if line in sample_taxa:
			s[line] = sample_taxa[line]
	
	try:
		primary,primary_counts = max(s.items(),key=lambda x: x[1][-1])
		tmp_primary,tmp_counts = max(s.items(),key=lambda x: x[1][-1])
		for z,w in s.items():
			if z != primary:
				shared,different = [0.0,0.0],[0.0,0.0] # [number_genes,number_reads]
				for j in range(255):
					if (mg_names[j] in MG_per_species[primary]) and (mg_names[j] in MG_per_species[z]):
						shared[0] += 1
						shared[1] += w[j]
					elif mg_names[j] in MG_per_species[z]:
						different[0] += 1
						different[1] += w[j]
				if len(MG_per_species[z]) < 5: # if taxa has les than 5 marker genes, remove
					del sample_taxa[z]
				elif 0.0 in shared:
					tmp_primary = tmp_primary.replace(' ','_')
					tmp_primary += '_???_'+z.split(';')[-1].replace(' ','_')
					tmp_counts += w
					del sample_taxa[z]
				elif (different[1]/shared[1]) > ((different[0]/shared[0])+0.1*different[0]/shared[0]): # a false positive would be expected to have a lower number of reads mapping to shared genes than to not-shared genes; here within an epsilon of 10% is allowed
					del sample_taxa[z]
					tmp_counts += w
				del sample_taxa[primary]
				sample_taxa[tmp_primary] = tmp_counts
	except: continue

aggregate_taxa = {}

for k,v in sample_taxa.items():
    if k == 'names': continue
    counts = v[-1] # number of mapped reads
    buscos = {mg_names[h] for h,i in enumerate(v) if i > 0} # number of mapped marker genes
    '''if counts < 2:
        if k in final_read_classifications: del final_read_classifications[k]
        continue
    if len(buscos) < 2:
        if k in final_read_classifications: del final_read_classifications[k]
        continue'''
    line = ''
    for j in k.split(';'):
        line += j+';'
        if line not in aggregate_taxa:
            aggregate_taxa[line] = [counts,buscos]
        else:
            aggregate_taxa[line][0] += counts
            aggregate_taxa[line][1].update(buscos)

with open('final_read_classifications.txt','w') as out:
	for k,v in final_read_classifications.items():
		for i in v:
			out.write(i+'\t'+k+'\n')

print('aggregated results. writing out to taxonomic report')

with open('Taxonomic_report.txt','w') as out:
	for k,v in sorted(aggregate_taxa.items()):
		counts,buscos = v[0],v[1]
		out.write(k+'\t'+str(counts)+'\t'+str(len(buscos))+'\n')





