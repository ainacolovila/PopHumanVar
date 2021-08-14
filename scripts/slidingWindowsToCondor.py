import pandas as pd
import numpy as np
import glob
from tqdm import tqdm
from multiprocessing import Pool
from plotnine import *

df = pd.read_csv("human_g1k_v37.fasta.fai",sep="\t",header=None,names=['chr','l','x1','x2','x3']).iloc[0:22,0:2]

w = 3*10**6

out = []
for index,row in df.iterrows():
	n = np.arange((row.l % w) /2,row.l + 1)
	tmp = np.lib.stride_tricks.sliding_window_view(n,3*10**6)[::10**6][:,[0,-1]].astype(int).astype(str)
	tmp = np.hstack([np.tile(row.chr,(tmp.shape[0],1)),tmp])

	toCondor = pd.DataFrame(tmp).apply(lambda row: row[0] + ":" + row[1] + "-" +row[2],axis=1)
	out.append(pd.DataFrame(tmp))

out = pd.concat(out)
run1 = out.iloc[0:1999,:]
run2 = out.iloc[1999:,:]

out.to_csv("/home/jmurga/Downloads/scripts/isafe_all_windows.txt",index=False,header=False,sep=',')
run1.to_csv("/home/jmurga/Downloads/scripts/isafe_windows_1.txt",index=False,header=False,sep=',')
run2.to_csv("/home/jmurga/Downloads/scripts/isafe_windows_2.txt",index=False,header=False,sep=',')