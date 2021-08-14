import pandas as pd
import numpy as np
from multiprocessing import Pool
import sys
import time 

def opening_isafe(path,population,df,nthreads):

	l = []
	for index, row in df.iterrows():
		try:
			tmp = pd.read_csv(path + "/" +population + "_" + str(row.chr) + "_" + str(row.start) + "_" + str(row.end) + ".iSAFE.out.gz",sep="\t")
			tmp['CHR'] = int(row.chr)
			l.append(tmp)
		except:
			next

	n = pd.concat(l).to_numpy()

	isafe = []
	for i in np.arange(21,23):
		start = time.time()
		print(i)
		tmp = n[n[:,3] == i]

		# only overlaping positions
		pos,u = np.unique(tmp[:,0],return_counts=True)
		s1 = pos[u > 2]

		pool = Pool(processes = nthreads)
		out1 = pool.starmap(genomewide_isafe,zip(s1,[tmp]*pos.shape[0]))
		pool.terminate()

		out1 = pd.DataFrame(np.vstack(out1),columns=['pos','isafe','daf','chr'])
		
		s2   = pd.DataFrame(pos[u < 3],columns=['pos'])
		daf  = pd.DataFrame({'pos':tmp[:,0],'isafe':np.nan,'daf':tmp[:,2],'chr':i}).drop_duplicates()
		out2 = pd.merge(s2,daf,how='left')
		out  = pd.concat([out1,out2]).sort_values('pos')

		isafe.append(out)
		print(time.time() - start)

	return(pd.concat(isafe))

def genomewide_isafe(p,m):
	flt = m[m[:,0] == p]
	return(flt[1])

path     = "/data/shared/isafe"
output   = "/home/jmurga/isafe"
population = sys.argv[1]

nthreads = int(sys.argv[2])

df       = pd.read_csv("/data/shared/isafe/isafe_all_windows.txt",header=None,names=['chr','start','end'])

tmp = opening_isafe(path + "/" + population,population,df,nthreads)

tmp = np.round(tmp,5)
tmpDf = pd.DataFrame(tmp,columns=['pos','isafe','daf','chr'])
tmpDf['population'] = population
tmpDf.to_csv(output + "/" + "isafe_" + population + ".txt.gz",compression="gzip")

########################def f(x):
import pandas as pd
import numpy as np
from multiprocessing import Pool
import sys
import time 

def opening_isafe(path,population,df,nthreads):

	l = []
	for index, row in df.iterrows():
		try:
			tmp = pd.read_csv(path + "/" +population + "_" + str(row.chr) + "_" + str(row.start) + "_" + str(row.end) + ".iSAFE.out.gz",sep="\t")
			tmp['CHR'] = int(row.chr)
			l.append(tmp)
		except:
			next

	n = pd.concat(l).to_numpy()

	isafe = []

	x = [pd.DataFrame(n[n[:,3] == i],columns=['pos','isafe','daf','chr']) for i in np.arange(1,23)]

	pool = Pool(processes = nthreads)
	out = pool.starmap(f,zip(x))
	pool.terminate()

	return(pd.concat(out).reset_index(drop=True))

def f(x):
	return x.groupby(['chr','pos']).isafe.apply(lambda r: r.values[1] if r.shape[0] > 2 else np.nan).reset_index()

path     = "/home/acolomer/getLDh/isafe"
output   = "/home/acolomer/getLDh/isafe"
# population = sys.argv[1]

nthreads = int(22)

df       = pd.read_csv("/home/acolomer/getLDh/isafe/isafe_all_windows.txt",header=None,names=['chr','start','end'])

population = np.array(["ACB","ASW","BEB","CDX","CEU","CHB","CHS","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","PJL","STU","TSI","YRI"])

for p in population:
	print(p)
	tmp = opening_isafe(path + "/" + p,p,df,nthreads)

	tmp = np.round(tmp,5)
	tmp['population'] = p
	tmp.to_csv(output + "/" + p + ".isafe.txt.gz",compression="gzip",na_rep="NA",index=False,sep="\t")
