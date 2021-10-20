import argparse
import numpy as np
from sklearn.linear_model import LogisticRegression

def train_logistic_regression(x):
      classifiers = {}
      for i in 'sfg':
          if i not in x: continue
          clf = LogisticRegression(C=1e5)
          v = np.array(x[i]+x['n'])
          vlabels = [1,1,0,0]
          V = v[:, np.newaxis]
          clf.fit(V, vlabels)
          classifiers[i] = [clf.coef_[0][0],clf.intercept_[0]]
      return classifiers

parser = argparse.ArgumentParser()
parser.add_argument("-s", help="strict classifier file")
args = parser.parse_args()

f = args.s   # 'strict_classifiers.txt'

with open(args.s+'.classifiers','w') as out:
	for h,i in enumerate(open(f)):
		line = i.strip().split('\t')
		prot = '_'.join(line[0].split('_')[:2])
		start = int(line[0].split('_')[-1])
		n,smin,smax,gmin,gmax,fmin,fmax = float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),float(line[7])
		if smin == 100: continue # skip regions with no species-specific read mappings
		if gmin == 100: gmin = 0
		if fmin == 100: fmin = 0
		if smax > max([n,gmin,gmax,fmin,fmax]):
			train = {}
			SMIN = max(smin,gmax,gmin,fmax,fmin,n)
			train['s'] = [smax,SMIN]
			if fmax > n:
				fmin = max(fmin,n)
				train['f'] = [fmax,fmin]
			if gmax > max(n,fmax):
				gmin = max(gmin,n,fmax)
				train['g'] = [gmax,gmin]
			non_negative_minimum = min(j for j in [fmin,gmin,SMIN] if j != 0)
			nmax = max(n,non_negative_minimum*0.8)
			nmin = nmax*0.8
			train['n'] = [nmax,nmin]
			results = train_logistic_regression(train)
			out.write(line[0]+'\t'+str(results)+'\n')


