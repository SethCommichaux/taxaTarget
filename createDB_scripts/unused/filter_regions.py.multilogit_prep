import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", help="strict classifier file")
args = parser.parse_args()

f = args.s   # 'strict_classifiers.txt'

classifiers = {}

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
			SMIN = max(smin,gmax,gmin,fmax,fmin,n)
			x = [str(smax),str(SMIN)]
			labels = ['s','s']
			if fmax > n:
				fmin = max(fmin,n)
				x += [str(fmax),str(fmin)]
				labels += ['f','f']
			if gmax > max(n,fmax):
				gmin = max(gmin,n,fmax)
				x += [str(gmax),str(gmin)]
				labels += ['g','g']
			non_negative_minimum = min(j for j in [fmin,gmin,SMIN] if j != 0)
			nmax = max(n,non_negative_minimum*0.8)
			nmin = nmax*0.8
			x += [str(nmax),str(nmin)]
			labels += ['n','n']
			dt = (n+non_negative_minimum)/2. # decision threshold
			out.write(line[0]+'\t'+','.join(x)+'\t'+','.join(labels)+'\t'+str(dt)+'\n')


