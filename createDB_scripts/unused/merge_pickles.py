import os
import pickle

results = {}

for i in os.listdir('.'):
	if i.endswith('pickles'):
		with open(i,'rb') as f:
			data = pickle.load(f)
		results.update(data)

with open('pickles','wb') as out:
	pickle.dump(results,out)
