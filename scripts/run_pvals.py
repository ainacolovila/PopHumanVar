import sys
import os
import pandas as pd 
import numpy as np
import glob
from multiprocessing import Pool
from tqdm import tqdm

# Which populations
path="/home/acolomer/phvDatabases"
test="ihs"
nthreads=20
files = glob.glob(path + "/" + test +"/*txt")

def open_df(x):
	return(pd.read_csv(x,header=0,sep='\t',na_values='\\N'))

pool = Pool(processes = nthreads)
out = pool.starmap(open_df,zip(files))
pool.terminate()

df = pd.concat(out)

pops = np.array(df.columns[3:])

chr_pos = df.iloc[:,0:3].values

out_pval = []
out_cutoff = []
pval_cutoff = 0.005

for p in tqdm(pops):

	tmp = pd.DataFrame({'chr':df.chr,'physicalPos':df.physicalPos,'rsid':df.rsid,p:df.loc[:,p]})
	stat = tmp.loc[:,p].values[~np.isnan(tmp.loc[:,p].values)]

	npCounts = np.vstack(np.unique(np.abs(stat),return_counts=True))
	denom = np.sum(npCounts[1])

	# Define function to get pvalue
	def pval(x,data=npCounts,d=denom):
		return(np.sum(data[1,data[0] >= x])/d)

	# Vectorize function so as it can be applyed through all ihs values
	vfunc = np.vectorize(pval)

	# Run the function
	pval = vfunc(npCounts[0])

	# Stack to ihs & pvalues and transpose
	stat_pval = np.vstack((npCounts[0], pval)).T

	# Get back to pandas, add suitable column names
	stat_pval_df = pd.DataFrame(data=stat_pval, columns=[p,'pvals'])

	df_pval_fill = pd.merge(tmp,stat_pval_df,how='left',on=[p])
 
	out_pval.append(df_pval_fill.pvals.values)

	out_cutoff.append(stat_pval[stat_pval[:,1] < pval_cutoff,0].min())

print("============> Stacking output")
df_pval = pd.DataFrame(np.hstack([chr_pos, np.vstack(out_pval).T]),columns = df.columns)
df_cutoff = pd.DataFrame(np.vstack(out_cutoff).T,columns = pops)

df_cutoff.insert(0, 'test', test)

print("============> Writting df's")
df_pval.to_csv(path + "/" + test + "/" + test + "_pvalues.txt.gz",sep='\t',index=False,na_rep=
'\\N')
df_cutoff.to_csv(path + "/" + test + "/" + test + "_cutoff.txt.gz",sep='\t',index=False)
(END)
