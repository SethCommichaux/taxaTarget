import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", help="strict classifier file")
parser.add_argument("-o", help="desired output file")
args = parser.parse_args()

f = args.s   # 'strict_classifiers.txt'

with open(args.o,'w') as out:
	for h,i in enumerate(open(f)):
		if h%10000 == 0: print(h)
		line = i.strip().split('\t')
		n,smin,smax,gmin,gmax,fmin,fmax = float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),float(line[7])
		if smin == 100: continue # skip regions with no species-specific read mappings
		if gmin == 100: gmin = 0
		if fmin == 100: fmin = 0
		if smax > max([n,gmax,fmax]):
			SMIN = max(smin,gmax,gmin,fmax,fmin,n)
			tmp = [n,SMIN,smax]
			if gmax > max(n,fmax):
				GMIN = max(gmin,n,fmax)
				tmp.append(GMIN)
				tmp.append(gmax)
			else:
				tmp += ['na','na']
			if fmax > n:
				FMIN = max(fmin,n)
				tmp.append(FMIN)
				tmp.append(fmax)
			else:
				tmp += ['na','na']
			tmp = [str(x) for x in tmp]
			out.write(line[0]+'\t'+'\t'.join(tmp)+'\n')	
