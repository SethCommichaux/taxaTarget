import argparse
import ast
import math

def logistic_function(x,x0,k,L=1):
	if x0 == 0: return 0
	prob = L/(1 + math.e**(-k*(x-x0)))	
	return prob

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-c", help="marker gene thresholds file")
parser.add_argument("-t", help="probability threshold",default=0.8,type=float)
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
	if classifiers[MG_region] != {}:
		mean_bitscore = bitscore/aln_len
		smin = classifiers[MG_region].get('s',0)
		gmin = classifiers[MG_region].get('g',0)
		fmin = classifiers[MG_region].get('f',0)
		nmax = classifiers[MG_region].get('n',0)	
		if nmax == 0:
			nmax = min([x for x in [smin,gmin,fmin] if x != 0])/2.0
		sprob = 0
		gprob = 0
		fprob = 0
		try:
			sprob = logistic_function(mean_bitscore,(smin+nmax)/2.0,k = abs(smin-nmax))
			if gmin != 0:
				gprob = logistic_function(mean_bitscore,(gmin+nmax)/2.0,k = abs(gmin-nmax))
			if fmin != 0:
				fprob = logistic_function(mean_bitscore,(fmin+nmax)/2.0,k = abs(fmin-nmax))
		except ZeroDivisionError: continue
		# nprob = logistic_function(mean_bitscore,nmax)
		# if nprob != 0: nprob = 1 - nprob
		probs = [fprob,gprob,sprob]
		ranks = 'fgs'
		rank = 'n'
		for x,prob in enumerate(probs):
			if prob != 0:
				if prob >= args.t:
					rank = ranks[x]
				else:
					break
		if rank != 'n':
			print(read,rank,probs,classifiers[MG_region],mean_bitscore)





'''


		elif mean_bitscore >= classifiers[MG_region]['s']: print(read,'s',mean_bitscore,classifiers[MG_region]['s'])
		elif mean_bitscore >= classifiers[MG_region].get('g',100): print(read,'g',mean_bitscore,classifiers[MG_region].get('g',100))
		elif mean_bitscore >= classifiers[MG_region].get('f',100): print(read,'f',mean_bitscore,classifiers[MG_region].get('f',100))
'''
