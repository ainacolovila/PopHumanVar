# a mano desde andromeda. Darwin no puede hacer este merge, guardar tabla y cargar en 
out = list()

for(i in 1:22){
	print(i)
	
	fnc = fread(paste0("/home/acolomer/phvDatabases/functional/functional_phv_",i,".txt"),na.strings='\\N')
	isafe = fread(paste0("/home/acolomer/phvDatabases/isafe/isafe_chr",i,".txt"),na.strings='\\N')
	nsl = fread(paste0("/home/acolomer/phvDatabases/nsl/nsl_chr",i,".txt"),na.strings='\\N')
	ihs = fread(paste0("/home/acolomer/phvDatabases/ihs/ihs_chr",i,".txt"),na.strings='\\N')
	geva = fread(paste0("/home/acolomer/phvDatabases/atlas_variant_age/atlas_chr",i,".txt.gz"),na.strings='\\N')

	######
	a = range(fnc$physicalPos) %>% t %>% as.data.table

	b = range(isafe$physicalPos) %>% t %>% as.data.table
	
	c = range(nsl$physicalPos) %>% t %>% as.data.table

	d = range(ihs$physicalPos) %>% t %>% as.data.table
	
	e= range(geva$physicalPos) %>% t %>% as.data.table

	x = rbind(a,b,c,d,e)

	chrom = data.table(chr=i,variable="chromRange",minV=min(x[,1]),maxV=max(x[,2]))

	########## functional
	f = fnc[,c('physicalPos','accession','NofPmids')]

	f$GWASAssocCounts <- NA

	if(sum(!is.na(f$accession))>0){
		f$NumGWASStudies[!is.na(f$accession)] <- do.call(rbind, lapply(strsplit(f$accession[!is.na(fnc$accession)],";"), function(x) length(x)))	
	}

	mn <- f[,c('NumGWASStudies','NofPmids')] %>% summarize_all(min,na.rm=T) 
	mx <- f[,c('NumGWASStudies','NofPmids')] %>% summarize_all(max,na.rm=T) 
	f = data.table(chr=rep(i,2),variable=names(mn),minV=as.numeric(mn),maxV=as.numeric(mx))


	###################
	geva_range = apply(geva[,c("AgeMode_Mut","AgeMode_Rec","AgeMode_Jnt")],2,range) %>% t %>% as.data.table

	tmp = data.table(chr=rep(i,3),variable=c("AgeMode_Mut","AgeMode_Rec","AgeMode_Jnt"),minV=geva_range$V1,maxV=geva_range$V2)

	###################

	out[[paste0(i)]] = rbind(chrom,f,tmp)
}

newslider = rbindlist(out)
newslider$variable = factor( as.character(newslider$variable), levels=c('GWASAssocCounts','maxNofPmids','AgeMode_Mut','AgeMode_Rec','AgeMode_Jnt','chromRange'))
newslider = newslider[order(variable,chr)]
fwrite(newslider,"/home/selector/sliderStaticValues.tsv",quote=F,col.names=F,sep='\t',na='\\N')

}
