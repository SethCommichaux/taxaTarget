import argparse
import ast
import math
from scipy.stats import norm

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-c", help="marker gene thresholds file")
args = parser.parse_args()

classifiers = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	classifiers[MG_region] = {}

for i in open(args.c):
	tmp = i.strip().split('\t')
	region = tmp[0]
	if region in classifiers:
		thresholds = ast.literal_eval(tmp[1])
		classifiers[region] = thresholds

for i in open(args.d):
	tmp = i.strip().split('\t')
	read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	probs=[0,0,0,0] # order = S,G,F,N
	if classifiers[MG_region] != {}:
		mean_bitscore = bitscore/aln_len
		if mean_bitscore <= classifiers[MG_region].get('n',0): continue
		else:
			smean = classifiers[MG_region]['s']
			sstd = 1.5 - (classifiers[MG_region]['s'] - max(classifiers[MG_region].get('g',0),classifiers[MG_region].get('f',0),classifiers[MG_region].get('n',0)))
			sprob = norm(smean,sstd).cdf(mean_bitscore)
			probs[0] = sprob
			if 'n' in classifiers[MG_region]:
				nmean = max(classifiers[MG_region]['n'],0.001)
				nstd = 1.5 - (min(classifiers[MG_region]['s'],classifiers[MG_region].get('g',100),classifiers[MG_region].get('f',100)) - classifiers[MG_region]['n'])
				nprob = 1 - norm(nmean,nstd).cdf(mean_bitscore)
				probs[3] = nprob
			if 'g' in classifiers[MG_region]:
				gmean = classifiers[MG_region]['g']
				gstd = 1.5 - (classifiers[MG_region]['s'] - classifiers[MG_region]['g'])
				gprob = norm(gmean,gstd).cdf(mean_bitscore)
				probs[1] = gprob
			if 'f' in classifiers[MG_region]:
				fmean = classifiers[MG_region]['f'] 
				fstd = 1.5 - (min(classifiers[MG_region]['s'],classifiers[MG_region].get('g',100)) - classifiers[MG_region]['f'])
				fprob = norm(fmean,fstd).cdf(mean_bitscore)
				probs[2] = fprob
			m = max(probs)
			if m >= 0.6:
				max_prob = ['s','g','f','n'][probs.index(m)]
				if max_prob != 'n':
					print(read,max_prob,probs,classifiers[MG_region],mean_bitscore)




'''


		elif mean_bitscore >= classifiers[MG_region]['s']: print(read,'s',mean_bitscore,classifiers[MG_region]['s'])
		elif mean_bitscore >= classifiers[MG_region].get('g',100): print(read,'g',mean_bitscore,classifiers[MG_region].get('g',100))
		elif mean_bitscore >= classifiers[MG_region].get('f',100): print(read,'f',mean_bitscore,classifiers[MG_region].get('f',100))
'''
