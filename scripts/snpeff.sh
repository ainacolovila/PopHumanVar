import pandas as pd
import numpy as np
import glob
from tqdm import tqdm
from multiprocessing import Pool

# def f(x,y):
#         z = x + ":" + y
#         return ",".join(np.unique(z))

# def u(x):
#         return ",".join(np.unique(x))

# out = []
# for n in np.arange(1,23):
	
# 	print(n)
	
# 	df = pd.read_csv("/data/shared/1000GP/annotations/snpeff/snpeff" + str(n) + ".txt",sep="\t",header=None,names=["CHROM","POS1","POS2","ALLELE","GENEID","FEATUREID","EFFECT","IMPACT"])

# 	tr = df.FEATUREID.apply(lambda row: np.char.array(row.split(','))).to_numpy()
# 	e  = df.EFFECT.apply(lambda row: np.char.array(row.split(','))).to_numpy()
# 	i  = df.IMPACT.apply(lambda row: np.char.array(row.split(','))).to_numpy()
# 	g  = df.GENEID.apply(lambda row: np.char.array(row.split(','))).to_numpy()

# 	pool = Pool(processes = 8)
# 	out1 = pool.starmap(f,zip(tr,e))
# 	pool.terminate()

# 	pool = Pool(processes = 8)
# 	out2 = pool.starmap(f,zip(g,e))
# 	pool.terminate()

# 	pool = Pool(processes = 8)
# 	out3 = pool.starmap(f,zip(g,i))
# 	pool.terminate()

# 	pool = Pool(processes = 8)
# 	out4 = pool.starmap(u,zip(e))
# 	pool.terminate()

# 	pool = Pool(processes = 8)
# 	out5 = pool.starmap(u,zip(i))
# 	pool.terminate()

# 	d = pd.DataFrame({'chr':df.CHROM,'pos':df.POS1,'gene:effect':out2,'transcript:effect':out1,'effect':out4,'gene:impact':out3,'impact':out5})
	
# 	d.to_csv("/data/shared/1000GP/annotations/snpeff/out"+str(n)+".txt",sep='\t',index=False)


###############

parallel -j10 -u "snpEff -Xmx8G ann -nostats -v -nodownload -geneId GRCh37.75 /data/shared/1000GP/chr{}_gp.vcf.gz | bcftools view - -Oz -o /data/shared/1000GP/annotations/snpeff/snpeff{}.vcf.gz" ::: `seq 1 22`

parallel -j8 -u "java -jar snpEff/SnpSift.jar extractFields snpeff{}.vcf.gz -s ',' -e '.' CHROM POS POS ID 'ANN[*].ALLELE' 'ANN[*].GENEID' 'ANN[*].FEATUREID' 'ANN[*].EFFECT' 'ANN[*].IMPACT' | intersectBed -a - -b /data/shared/1000GP/masks/chr{}_mask.bed > snpeff{}.txt" ::: `seq 1 22`


