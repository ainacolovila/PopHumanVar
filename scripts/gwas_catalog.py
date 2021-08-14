import pandas as pd
import numpy as np
from tqdm import tqdm

df = pd.read_csv("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",sep='\t')

df.columns = ["DATE_ADDED_TO_CATALOG","PUBMEDID","FIRST_AUTHOR","DATE","JOURNAL","LINK","STUDY","DISEASE_TRAIT","INITIAL_SAMPLE_SIZE","REPLICATION_SAMPLE_SIZE", "REGION","CHR_ID","CHR_POS","REPORTED_GENEs","MAPPED_GENE","UPSTREAM_GENE_ID","DOWNSTREAM_GENE_ID","SNP_GENE_IDS","UPSTREAM_GENE_DISTANCE","DOWNSTREAM_GENE_DISTANCE","STRONGEST_SNP_RISK_ALLELE","SNPS","MERGED","SNP_ID_CURRENT","CONTEXT","INTERGENIC","RISK_ALLELE_FREQUENCY","PVALUE","PVALUE_MLOG","PVALUE_text","OR_or_BETA","95_CI","PLATFORM_SNPS_PASSING_QC","CNV","MAPPED_TRAIT","MAPPED_TRAIT_URI","STUDY_ACCESSION","GENOTYPING_TECHNOLOGY"]

tmp = df[['CHR_ID','CHR_POS','STRONGEST_SNP_RISK_ALLELE','STRONGEST_SNP_RISK_ALLELE','RISK_ALLELE_FREQUENCY','DISEASE_TRAIT','STUDY_ACCESSION','CONTEXT']]

tmp.columns = ["chr","physicalPos",'rsid','risk_allele',"risk_allele_freq","trait","accession","context"]

tmp.chr = tmp.chr.astype(str)
tmp.trait = tmp.trait.astype(str)
tmp.risk_allele_freq = tmp.risk_allele_freq.astype(str)

tmp = tmp[tmp.chr != 'nan']
tmp = tmp[~tmp.chr.str.contains('x')]
tmp = tmp[(tmp.chr !='X') & (tmp.chr != 'Y')]
tmp = tmp.reset_index(drop=True)

def flt(x,n):
    try: 
        if ';' in x:
            return x.split(';')
        else:
            return [x] * n 
    except: 
        return [x] * n 

out = list()
for index,row in tqdm(tmp.iterrows(),total=tmp.shape[0]):
    if(';' in row.chr):
        n = row.chr.count(';') + 1
        x = row.apply(lambda row: flt(row,n))
        for z in np.arange(0,n):
            out.append(pd.DataFrame([[x.chr[z],x.physicalPos[z],x.rsid[z],x.risk_allele[z],x.risk_allele_freq[z],x.trait[z],x.accession[z],x.context[z]]],columns=x.index.values))
    else:
        out.append(tmp.iloc[[index],:])

gwas_catalog = pd.concat(out)
gwas_catalog.chr = gwas_catalog.chr.astype(int)
gwas_catalog.physicalPos = gwas_catalog.physicalPos.astype(int)
gwas_catalog.rsid = gwas_catalog.rsid.str.split('-', 1, expand=True).iloc[:,0]
gwas_catalog.risk_allele = gwas_catalog.risk_allele.str.split('-', 1, expand=True).iloc[:,1]
gwas_catalog.risk_allele_freq[gwas_catalog.risk_allele_freq == 'NR'] = 'NA'
gwas_catalog.risk_allele_freq[gwas_catalog.risk_allele_freq == 'nan'] = 'NA'

gwas_catalog.to_csv("/home/acolomer/phvDatabases/gwas_catalog/gwas_to_r.txt",index=False,header=True,sep='\t',na_rep='NA')

#######
df = fread("gwas_to_r.txt")
gwas_catalog_collapse = df %>% group_by(chr,physicalPos,rsid) %>% summarize("risk_allele"=paste0(unique(risk_allele),collapse=';'), "risk_allele_freq"=paste0(unique(risk_allele_freq),collapse=';'), "trait"=paste0(unique(trait),collapse=';'),"accession"=paste0(unique(accession),collapse=';'),"context"=paste0(unique(context),collapse=';')) %>% as.data.table

gwas_pos = gwas_catalog_collapse[,1:3]
gwas_pos[['end']] = gwas_pos[['physicalPos']]
gwas_pos[['physicalPos']] = gwas_pos[['physicalPos']] - 1
gwas_pos = gwas_pos[,c(1,2,4,3)]
gwas_pos[['chr']] = paste0('chr', gwas_pos[['chr']])

fwrite(gwas_pos,file="gwas_catalog_positions_hg38.txt",sep='\t',col.names=F)

#liftOver /home/acolomer/phvDatabases/gwas_catalog/gwas_catalog_positions_hg38.txt /home/acolomer/phvDatabases/gwas_catalog/hg38ToHg19.over.chain.gz /home/acolomer/phvDatabases/gwas_catalog/gwas_catalog_positions_hg19.txt unMapped.txt && cut -f1,2,3,4 gwas_catalog_positions_hg19.txt | sed 's/chr//g' > tmp && mv tmp  gwas_catalog_positions_hg19.txt


hg38 = fread("sed 's/chr//g' gwas_catalog_positions_hg38.txt",col.names=c("chr",'startOld','physicalPos','rsid'))
hg19 = fread("gwas_catalog_positions_hg19.txt",col.names=c("chr",'start','end','rsid'))

liftOver = merge(hg19[,c(1,4,3)],hg38)[,c(1,5,2,3)]

x = merge(gwas_catalog_collapse,liftOver,by=c("chr","physicalPos","rsid"))
x = x[,c("chr","end",'rsid','risk_allele',"risk_allele_freq","trait","accession","context")]
names(x)[2] = "physicalPos"
x = x[order(chr,physicalPos)]

fwrite(x,file="gwas_catalog.txt",sep='\t',na='\\N',quote=F)

# BEFORE PUSH
# sed -i '1d' gwas_catalog.txt 

# sed 's/NA;//g' gwas_catalog.txt  | sed 's/;NA//g' > tmp && mv tmp gwas_catalog.txt 