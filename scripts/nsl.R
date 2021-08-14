library(data.table)
library(fst)

out = list()
output = ""
# nsl too
for(f in list.files("/data/shared/ihs/",pattern="norm.gz",full.names=T)){
	print(f)
	tmp = fread(f)

	name = unlist(strsplit(f,'/'))
	name = name[length(name)]
	population = unlist(strsplit(name,"\\."))

	tmp[['population']] = population[1]

	out[[paste0(f)]] = tmp[,c(1,2,3,8,10)]
}

df = rbindlist(out)

for(i in 1:22){
	print(i)
	x = dcast.data.table(df[chr==i],chr + physicalPos + rsid  ~ population,value.var = "standarizediHS")
	fwrite(x,file=paste0(output,"/ihs_chr",i,".txt"),sep='\t',na='\\N',row.names=F,col.names=T,quote=F)
}

####################NSL

populations = c("ACB","ASW","BEB","CDX","CEU","CHB","CHS","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","PJL","STU","TSI","YRI")
path = "/home/acolomer/getLDh/nSL/normXPEHH/"

for(p in populations){
	out = list()
	print(p)
	for(i in 1:22){
		print(i)
		tmp = fread(paste0(path,p,"/chr",i,"_",p,".nsl.out.100bins.norm"))

		names(tmp) = c("rsid","physicalPos","1freq","sl1","sl0","unstandarizednSL","standarizednSL","extremevalue")
		tmp = cbind(chr=i,tmp)
		out[[paste0(i)]] = tmp
	}
	df = rbindlist(out)
	fwrite(df,file=paste0(path,"/",p,".nsl.norm.gz"),sep='\t',row.names=F)
}

####################ISAFE

path = "/data/shared/isafe/tmp"
output = "/data/shared/isafe/"
	
rs = fread("/home/acolomer/getLDh/Data/SNPs/rsid_bialleics_ancestral.txt.gz")
names(rs) = c('chr','physicalPos','rsid')
# rs[['chr']] = as.numeric(rs[['chr']])


for(i in list.files(path,pattern='gz',full.names=T)){

	print(i)

	tmp = fread(i)
	names(tmp)[2] = "physicalPos"
	df  = merge(tmp,rs,all.x=T)

	df = df[,c(1,2,5,3,4)]
	x = unlist(strsplit(i,"/"))
	name = x[length(x)]

	fwrite(df,file=paste0(output,"/",name),row.names=F,sep='\t',na="NA",quote=F)

}

library(data.table)
library(fst)

out = list()
output ="/home/acolomer/phvDatabases/isafe"

for(f in list.files("/data/shared/isafe/",pattern="txt.gz",full.names=T)){
	print(f)
	tmp = fread(f)

	out[[paste0(f)]] = tmp
}

df = rbindlist(out)

for(i in 1:22){
	print(i)
	x = dcast.data.table(df[chr==i],chr + physicalPos + rsid  ~ population,value.var = "isafe")

	idx = rowSums(x[,4:25],na.rm=T) != 0

	fwrite(x,file=paste0(output,"/isafe_chr",i,".txt"),sep='\t',na='\\N',row.names=F,col.names=T,quote=F)
}


for(i in list.files()){
	print(i)
	x = fread(i,na.strings='\\N')
	idx = rowSums(x[,4:25],na.rm=T) != 0
	y = x[idx,]
	fwrite(y,file=paste0("/home/acolomer/phvDatabases/isafe/",i),sep='\t',na='\\N',row.names=F,col.names=T,quote=F)

}

###################GEVA
library (tidyr)
library(dplyr)
library(data.table)


for (n in 1:22) {
  print(paste0('-->', n))
   
   Age <- fread(file = paste0('/home/acolomer/phvDatabases/atlas_variant_age/rawData/atlas.chr',n,'.csv.gz'), header = TRUE, sep = ',')

  
  colnames(Age)[1:3] <- c('rsid','chr','physicalPos')
  tmp = cbind(Age[,c(3,1)],Age[,4:27])
  Age_fltr <- tmp[DataSource == 'TGP']  %>% select(physicalPos, rsid, AlleleRef, AlleleAlt, AlleleAnc, starts_with('Age'), starts_with('Qual')) %>% as.data.table
  
  fwrite(Age_fltr, file = paste0('/home/acolomer/phvDatabases/atlas_variant_age/atlas_chr',n,'.txt'),sep = "\t", na='\\N',quote=F)
 
}

###################SNPEFF

for(i in 1:22){
	print(i)
	df = fread(paste0("/data/shared/1000GP/annotations/snpeff/snpeff",i,".txt.gz"))
	tmp = df[,c(1,3,4,6,8,9)]
	names(tmp) = c('chr','physicalPos','rsid','gene','effect','impact')
	tmp = tmp[grepl("^rs", tmp$rsid)]

	fwrite(tmp,file=paste0("/home/acolomer/phvDatabases/snpeff/snpeff_phv_",i,".txt.gz"),na='\\N',quote=F,sep='\t')

}

###################DISGENET
#Hg38To19
df = fread("/home/acolomer/phvDatabases/disgenet/rawData/all_variant_disease_associations.tsv.gz")

tmp = df[chromosome!='Y' & chromosome!='X' & chromosome!='MT']

z = tmp %>% group_by(snpId,chromosome,position) %>% summarize("diseaseName" = paste0(unique(diseaseName),collapse=';'),"diseaseId" = paste0(unique(diseaseId),collapse=';'),"diseaseType" = paste0(unique(diseaseType),collapse=';'),"source" = paste0(unique(source),collapse=';'),"DSI" = paste0(unique(DSI),collapse=';'),"DPI" = paste0(unique(DPI),collapse=';'),"EI" = paste0(unique(EI),collapse=';'),"NofPmids" = paste0(unique(NofPmids),collapse=';')) %>% as.data.table

x = z[,2:3]
names(x) = c("chr","start")

x$end = x$start
x$start = x$start - 1
x$chr = paste0('chr',x$chr)
x$id = z[,1]

#bedfile
fwrite(x,file="/home/acolomer/phvDatabases/disgenet/rawData/disgenet_positions_hg38.txt",col.names=F,sep='\t')

# liftOver /home/acolomer/phvDatabases/disgenet/rawData/disgenet_positions_hg38.txt /home/acolomer/phvDatabases/disgenet/rawData/hg38ToHg19.over.chain.gz /home/acolomer/phvDatabases/disgenet/rawData/disgenet_positions_hg19.txt unMapped.txt

disgenet19 = fread("/home/acolomer/phvDatabases/disgenet/rawData/disgenet_positions_hg19.txt",col.names=c('chromosome','start','end','snpId'))
disgenet19$chromosome = gsub("chr","",disgenet19$chromosome)

y = left_join(z,disgenet19,by=c("chromosome","snpId")) %>% as.data.table
y = merge(z,disgenet19,by=c("chromosome","snpId")) %>% as.data.table
y = y[,c("chromosome","end","snpId","diseaseName","diseaseId","diseaseType","source","DSI","DPI","EI","NofPmids")]
names(y)[1:3] = c("chr","physicalPos","rsid")
y$physicalPos = as.numeric(y$physicalPos)
y$chr = as.numeric(y$chr)
y$DSI      = as.numeric(y$DSI)
y$DtxtPI      = as.numeric(y$DPI)
y$EI       = as.numeric(y$EI)
y$NofPmids = as.numeric(y$NofPmids)


y = y[order(chr,physicalPos)]
fwrite(y,file="/home/acolomer/phvDatabases/disgenet/disgenet.txt.gz",col.names=T,sep='\t',na='\\N',quote=F)

################################clinvar

df = fread("/home/acolomer/phvDatabases/clinvar/Final_Aina.txt")

tmp = df %>% group_by(chr,physicalPos,rsid) %>% summarize("ClinicalSignificance"=paste0(unique(ClinicalSignificance),collapse=';'),"ClinSigSimple"=paste0(unique(ClinSigSimple),collapse=';'),"PhenotypeList"=paste0(unique(PhenotypeList),collapse=';')) %>% as.data.table
fwrite(tmp,file="/home/acolomer/phvDatabases/clinvar/Final_Aina_Collapse.txt.gz",col.names=T,sep='\t',na='\\N',quote=F)

################################clinvar
df = fread("/home/acolomer/phvDatabases/regulomedb/rawData/TSTFF344324.tsv.gz")
df$chrom = as.numeric(gsub("chr","",df$chrom))

for(i in 1:22){

	tmp = df[chrom == i]
	tmp = tmp[,!"start"]
	names(tmp) = c('chr','physicalPos','rsid','TFbinding','DNasePeak','motif','DNaseFootprint','eQTL','matchedTFmotif','matchedDNaseFootprint','ranking')

	fwrite(tmp,file=paste0("/home/acolomer/phvDatabases/regulomedb/ReguDB_chr,",i,".txt"),col.names=T,sep='\t',na='\\N',quote=F)

}
################################merge all

g = fread("/home/acolomer/phvDatabases/gwas_catalog/gwas_catalog.txt.gz",na.strings='NA')
d = fread("/home/acolomer/phvDatabases/disgenet/disgenet.txt.gz",na.strings='\\N')
c = fread("/home/acolomer/phvDatabases/clinvar/Final_Aina_Collapse_Clinvar.txt.gz",na.strings='\\N')


for(i in 1:22){

	print(i)

	#snpeff
	s = fread(paste0('/home/acolomer/phvDatabases/snpeff/snpeff_phv_',i,'.txt.gz'),na.string='\\N')
	
	#regulome db
	r = fread(paste0("/home/acolomer/phvDatabases/regulomedb/ReguDB_chr",i,".txt.gz"))
	#gwas
	g1 = g[chr==i]
	
	#disgenet
	d1 = d[chr==i]

	#clinvar
	c1 = c[chr==i]

	x = Reduce(function(x, y) merge(x, y, all=TRUE), list(s,r,g1,d1,c1))
	
	# z = Reduce(function(x, y) merge(x, y), list(tmp,r,g1))
	
	# out1 = merge(tmp,r)
	# out2 = merge(out1,g1)
	x = x[,2:ncol(x)]
	x$TFbinding = as.numeric(x$TFbinding)
	x$DNasePeak = as.numeric(x$DNasePeak)
	x$motif = as.numeric(x$motif)
	x$DNaseFootprint = as.numeric(x$DNaseFootprint)
	x$eQTL = as.numeric(x$eQTL)
	x$matchedTFmotif = as.numeric(x$matchedTFmotif)
	x$matchedDNaseFootprint = as.numeric(x$matchedDNaseFootprint)

	fwrite(x,file=paste0('/home/acolomer/phvDatabases/functional/functional_phv_',i,'.txt'),na='\\N',sep='\t',quote=F,row.names=F)
}

###########
twoSplitDf = function(dfFiltered,c1=c("gene","impact","effect"),c2=c("effect"),s1=",",s2="&"){
    z1 = cSplit(dfFiltered,c1,sep=',',direction="long",makeEqual=F)
    if(!is.null(c2) & !is.null(s2)){
        z2 = cSplit(z1,c2,sep='&',direction="long") 
        return(z2)
    }else{
        return(z1)   
    }
}

##################CHROMRANGES
out = list()

for(i in 1:22){
	print(i)
	fnc = range(fread(paste0("/home/acolomer/phvDatabases/functional/functional_phv_",i,".txt"))$physicalPos) %>% t %>% as.data.table

	isafe = range(fread(paste0("/home/acolomer/phvDatabases/isafe/isafe_chr",i,".txt"))$physicalPos) %>% t %>% as.data.table
	
	nsl = range(fread(paste0("/home/acolomer/phvDatabases/nsl/nsl_chr",i,".txt"))$physicalPos) %>% t %>% as.data.table

	ihs = range(fread(paste0("/home/acolomer/phvDatabases/ihs/ihs_chr",i,".txt"))$physicalPos) %>% t %>% as.data.table

	x = rbind(fnc,isafe,nsl,ihs)

	y = data.table(chr=i,variable="chromRange",minV=min(x[,1]),maxV=max(x[,2]))
	out[[paste0(i)]] = y
}


df = rbindlist(out)

#######################PHS
df = fread("https://pophumanscan.uab.cat/data/files/rawPophumanscanTable.tab")
df = df[chr!='chrX' & chr!='chrY']
df$chr = as.numeric(gsub("chr","",df$chr))
for(i in 1:22){

	print(i)
	tmp = df[chr==i]
	tmp = tmp[,c('start','end','GeneID','pops','source')]
	tmp = tmp[order(start,end)]
	fwrite(tmp,paste0("/home/acolomer/phvDatabases/phs/phs_chr",i,".txt"),sep='\t',na='\\N',quote=F)

}

###########################SLIDER
con <- dbConnect(RMySQL::MySQL(), user='shiny',password='***',dbname='functional',host='localhost')
	myQuery <- paste0("SELECT * FROM sliderStaticValues")
slider =  dbGetQuery(con, myQuery) %>% as.data.table
dbDisconnect(con)

# Connect & get SQL informaiton
out = list()
for(i in 1:22){
	print(i)
	slider_chr = slider[chr == i][c(1,2,6)]
	geva = sqlToDf(i,1, 10^9, 'geva') %>% select(starts_with('AgeMode'))
	names(geva)
	geva_range = apply(geva,2,range) %>% t %>% as.data.table
	tmp = data.table(chr=rep(i,3),variable=names(geva),minV=geva_range$V1,maxV=geva_range$V2)
	out[[paste0(i)]] = rbind(slider_chr,tmp)
}

newslider = rbindlist(out)
newslider$variable = factor( as.character(newslider$variable), levels=c('GWASAssocCounts','maxNofPmids','AgeMode_Mut','AgeMode_Rec','AgeMode_Jnt','chromRange'))
newslider = newslider[order(variable,chr)]
fwrite(newslider,"/home/selector/sliderStaticValues.tsv",quote=F,col.names=F,sep='\t',na='\\N')
#################STATQUANTILE

sqlToDf <- function(chr, coordStart, coordEnd, dbName){
        con <- dbConnect(RMySQL::MySQL(),user='shiny',password='***',dbname=dbName,host='localhost')
        myQuery <- paste0("SELECT * FROM ", "chr", chr, " WHERE physicalPos BETWEEN ", coordStart, " AND ", coordEnd)
        df <- dbGetQuery(con, myQuery)
        dbDisconnect(con)
        return(df %>% as.data.table)
}

collapseDB = function(stat,threshold=0.9999){

	out = list()

	for(i in 1:22){

		print(i)

		tmp1 = sqlToDf(i,1,10^20,stat)
		out[[i]] = tmp1

	}
	return(rbindlist(out))
}

ihs = collapseDB('ihs')
nsl = collapseDB('nsl')
isafe = collapseDB('isafe')

ihs = as.data.table(read.fst('ihs.fst'))
nsl = as.data.table(read.fst('nsl.fst'))
isafe = as.data.table(read.fst('isafe.fst'))

stat_quantile = function(df,pops=c("YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW", "CEU", "TSI", "FIN", "GBR", "IBS","CHB", "JPT", "CHS", "CDX", "KHV", "GIH", "PJL", "BEB", "STU", "ITU")
,threshold=0.9999){

	cutoff = list()
	for(p in pops){
		print(p)
		df_pop = df[[p]]
		df_pos = quantile(abs(na.omit(df_pop)),probs = c(1-(threshold/100))) %>% as.vector()
		# df_pos = quantile(na.omit(df_pop[df_pop>0]),probs = c(threshold)) %>% as.vector()
		# df_neg = quantile(na.omit(df_pop[df_pop < 0 ]),probs = c(1-threshold)) %>% as.vector()

		cutoff[[p]] = data.table(pop = p,cutoff=df_pos)
	}
	cutoff     = rbindlist(cutoff)
	return(cutoff)
}

cutoff_ihs = stat_quantile(df=ihs,threshold=0.9995)
write.fst(cutoff_ihs,"/home/selector/cutoff_ihs.fst")

cutoff_nsl = stat_quantile(df=nsl,threshold=0.9995)
write.fst(cutoff_nsl,"/home/selector/cutoff_nsl.fst")

cutoff_isafe = stat_quantile(df=isafe,pops=c('CEU'),threshold=0.01)
write.fst(cutoff_isafe,"/home/selector/cutoff_isafe.fst")