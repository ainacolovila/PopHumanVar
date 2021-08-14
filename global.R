## global.R ##

## app.R ##
# PACKAGES
suppressWarnings(library(shiny))
suppressWarnings(library(shinydashboard))
suppressWarnings(library(shinyBS))
suppressWarnings(library(shinyjs))
suppressWarnings(library(shinyalert))
suppressWarnings(library(shinyWidgets))
suppressWarnings(library(shinythemes))
suppressWarnings(library(shinydashboardPlus))
suppressWarnings(library(shinyhelper))
suppressWarnings(library(shinycustomloader))
suppressWarnings(library(shinycssloaders))
suppressWarnings(library(ggplot2))
suppressWarnings(library(plotly))
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyverse))
suppressWarnings(library(tibble))
suppressWarnings(library(data.table)) 
suppressWarnings(library(DBI))
suppressWarnings(library(dbplyr))
suppressWarnings(library(RMySQL))
suppressWarnings(library(splitstackshape))  
suppressWarnings(library(stringr))  
suppressWarnings(library(DT))
suppressWarnings(library(JBrowseR))

##################### FUNCTIONS

# Connect & get SQL informaiton
sqlToDf <- function(chr, coordStart, coordEnd, dbName){
	con <- dbConnect(RMySQL::MySQL(),
					 user='shiny',
					 password='***',
					 dbname=dbName,
					 host='localhost')
	myQuery <- paste0("SELECT * FROM ", 
					  "chr", chr, 
					  " WHERE physicalPos BETWEEN ", coordStart, 
					  " AND ", coordEnd)
	df <- dbGetQuery(con, myQuery)
	dbDisconnect(con)
	return(df %>% as.data.table)
}

# Connect phs and overalp

phs_overlap <- function(chr,coordStart,coordEnd){

	con <- dbConnect(RMySQL::MySQL(),
					 user='shiny',
					 password='***',
					 dbname='phs',
					 host='localhost')
	myQuery <- paste0("SELECT * FROM ", 
					  "chr", chr)
	phs <- as.data.table(dbGetQuery(con, myQuery))
	dbDisconnect(con)

	phs$chr = chr
	phs = phs[,c('chr','start','end')]
	setkey(phs,chr,start,end)

	tmp = data.table(chr=chr,start=coordStart,end=coordEnd)
	setkey(tmp,chr,start,end)
	
	phs_phv_overlap = foverlaps(tmp, phs, type="any", nomatch=0L)

	if(nrow(phs_phv_overlap) < 1){
		return(NULL)
	}else{
		s = c(phs_phv_overlap$start,phs_phv_overlap$i.start)
		e = c(phs_phv_overlap$end,phs_phv_overlap$i.end)
		return(paste0('chr',phs_phv_overlap$chr,':',min(s),'-',max(e)))
	}
}

# To get rsID coordinates
getID <- function(idInput){
	hg19 <- dbConnect(MySQL(), 
				   user="genome", 
				   db="hg19", 
				   host="genome-mysql.cse.ucsc.edu")
	myQuery <- paste0('SELECT chrom,chromStart,chromEnd FROM snp147 WHERE name = "',idInput, '"')
	coord <- dbGetQuery(hg19, myQuery) %>% suppressWarnings()
	dbDisconnect(hg19)
	
	names(coord) <- c("chr","start","end")
	coord$chr <- gsub("chr","",coord$chr)
	return(coord)
}

# Get coordinates from EnsemblID
getEnsID <- function(idInput){
	# Connection
	con <- dbConnect(RMySQL::MySQL(), 
				  user='shiny', 
				  password='***', 
				  dbname='functional', 
				  host='localhost')
	
	myQuery <- paste0('SELECT chr,start,end FROM genes WHERE ensemblID = "',idInput, '"')
	coord <- dbGetQuery(con, myQuery) %>% suppressWarnings()
	dbDisconnect(con)
	
	return(coord)
}

# Get coordinates from hugoGeneID
getHUGO <- function(idInput){
	# Connection
	con <- dbConnect(RMySQL::MySQL(), 
				  user='shiny', 
				  password='***', 
				  dbname='functional', 
				  host='localhost')

		myQuery <- paste0('SELECT chr,start,end FROM genes WHERE hgnc = "',idInput, '"')
	coord <- dbGetQuery(con, myQuery) %>% suppressWarnings()
	dbDisconnect(con)
	
	return(coord)
}

# Update max an min from widgets when chr change
getSliderValues <- function(chrN, clock){
	
	con <- dbConnect(RMySQL::MySQL(), 
		user='shiny',
		password='***',
		dbname='functional',host='localhost')
	myQuery <- paste0("SELECT * FROM sliderStaticValues")
	queryInfo <- dbGetQuery(con, myQuery) %>% as.data.table
	dbDisconnect(con)
	
	if(!is.null(clock)){
		# Modify data
		infoSlider <- queryInfo %>% dplyr::filter(chr %in% chrN, variable %in% c("GWASAssocCounts","maxNofPmids",paste0("AgeMode_",clock))) %>% dplyr::select(variable,minV,maxV)
	} else {
			infoSlider <- queryInfo %>% dplyr::filter(chr %in% chrN, variable=="chromRange") %>% 
			dplyr::select(minV,maxV)
	}
	return(infoSlider)
}

# Funció colorByQuality
colorByQuality <- function(min, max){
	colorScale <- c("saddlebrown","darkred","red","chocolate1","orange",
				 "goldenrod1", "yellow2","lightgreen", 
				 "green","forestgreen")
	calidade <- colorRampPalette(colorScale[min:max])
	return(calidade(10))
}

# Function high plots
returnHeigh <- function(selectedPops){
	return(sizesTab <- length(selectedPops)*200)
}

# Function high plots
returnHeighPop <- function(selectedPops){
	nPops <- length(selectedPops)
	if(nPops == 1){
		return(1)
	}
	else if(nPops<=6){
		lateral = 1/nPops
		center  = lateral
		t = nPops-2
		return(c(lateral,rep(center, t),lateral))
	} else {
    propStd = 1/nPops
    lateral = propStd - (propStd*0.28)
    center  = (1-(lateral*2))/(nPops-2)
    t = nPops-2
  	return(c(lateral,rep(center, t),lateral))
	}
}

# Function to search ID
getGenomicCoordinates <- function(idInput){
	if(grepl('^rs[[:digit:]]+', idInput)){
		return(getID(idInput))
	}
	if(grepl('^ENSG[[:digit:]]+', idInput)){
		return(getEnsID(idInput))
	}
	else{
		return(getHUGO(idInput))
	}
}

# Function to split collapsed dataframe
SplitDf <- function(dfFiltered, c1=c("gene_effect","impact","effect"),c2=c("effect"),s1=",",s2="&"){
	
	z1 = cSplit(dfFiltered, c1, sep=',',direction="long", makeEqual=F)

	if(!is.null(c2) & !is.null(s2)){
		z2 = cSplit(z1,c2,sep='&',direction="long")

		ind = which(names(z2) %in% c1,T)

		z2[, (ind) := lapply(.SD, as.character), .SDcols = ind]

		return(z2)

	}else{
		return(z1)  
	}
} 

# Merge all df for the region
mergeAllDF <- function(chr, coordStart, coordEnd){
	
	# Load data
	ihsDF <- sqlToDf(chr, coordStart, coordEnd, 'ihs')
	nslDF <- sqlToDf(chr, coordStart, coordEnd, 'nsl')
	functionalDF <- sqlToDf(chr, coordStart, coordEnd, 'functional')
	gevaDF <- sqlToDf(chr, coordStart, coordEnd, 'geva')
	isafeDF <- sqlToDf(chr, coordStart, coordEnd, 'isafe')

	# Merge All
	allDF = Reduce(function(x, y) merge(x, y, 
		all=TRUE, by=c("physicalPos", "rsid"), 
		suffixes = c("_x","_y")), 
	list(ihsDF, nslDF, functionalDF, gevaDF, isafeDF))
	
	# Change names
	names(allDF)<- gsub("_x", "_ihs", names(allDF))
	names(allDF)<- gsub("_y", "_nsl", names(allDF))
	names(allDF)[95:116] <- paste0(names(allDF)[95:116], "_isafe")
	
	allDF$NumGWASStudies <- NA
	allDF$NumDGNAsso <- NA

	if(sum(!is.na(allDF$accession))>0){
		allDF$NumGWASStudies[!is.na(allDF$accession)] <- do.call(rbind, lapply(strsplit(allDF$accession[!is.na(allDF$accession)],";"), function(x) length(x)))	
	}
	if(sum(!is.na(allDF$diseaseId))>0){
		allDF$NumDGNAsso[!is.na(allDF$diseaseId)] <- do.call(rbind, lapply(strsplit(allDF$diseaseId[!is.na(allDF$diseaseId)],";"), function(x) length(x)))
	}
	return(allDF)
}

# To melt the DF (plots separated by pop or metapop)
meltDFtoPlot <- function(df, varCol){
	if(varCol=='isafe'){
    prefiltered    = melt(setDT(df), id.vars = c("rsid", "physicalPos"), variable.name = "POP")
    melted         = prefiltered %>% filter(!is.na(value))
    melted$METAPOP = NA
		for(MM in metapopsData$meta){melted$METAPOP[melted$POP %in% get(paste0(MM,'pops'))]<- MM}
	}else{
			# filter <- df 
			prefiltered <- melt(setDT(df), id.vars = c("rsid", "physicalPos"), variable.name = "POP")
			# prefiltered$POP<- gsub(paste0("_",varCol),"",prefiltered$POP)
			melted <- prefiltered %>% filter(!is.na(value))
	}
	return(as.data.table(melted))
}

# To generate general boxplots
plotBoxPlot <- function(df_all, df_filtered, chrN, extremeValuesInput=FALSE, selectedPops,stat){
	df_filtered = df_filtered[thres_ad != "Not extreme"]
	fig <- df_filtered %>% 
		plot_ly(
			type = 'scatter',
			mode = 'markers', 
			x = ~POP,
			y = ~value, 
			color = ~POP, 
			colors= as.vector(popsData$colors[match(selectedPops,popsData$pops)]),
			text = ~paste0('<b>',rsid,'</b><br>chr', chrN,':', physicalPos),
			hovertemplate = paste(
				"%{text}<br>",
				"<b>",stat,"</b>: %{y}<br>",
				"<b>Population</b>: %{x}"
			),
			showlegend = FALSE)
	# if Filter by extreme, change boxplot colors
	if(extremeValuesInput==TRUE){
		fig <- fig %>% add_trace(
			data = df_all,
			y = ~value,
			type="box",
			boxpoints=FALSE,
			showlegend = FALSE,
			fillcolor ="#e7efef",
			line = list(color = "#909495"))
	}else{
		fig <- fig %>% add_trace(
			data = df_all,
			y = ~value,
			type="box",
			boxpoints=FALSE,
			showlegend = FALSE,
			color = ~POP,
			colors= as.vector(popsData$colors[match(selectedPops,popsData$pops)]))
	}
	fig <- fig %>% layout(xaxis = list(categoryorder = "array", categoryarray = as.vector(popsData$pops[popsData$pops %in% selectedPops])))
	return(fig)
}

# To generate bypop Plots
plotScatterPlot <- function(df, varStat,column_color,popVector){
	
	if(column_color == 'AD'){
		pal = c('firebrick', 'darkblue')
	}else{
		pal <- c('firebrick','darkblue','gray')
		pal <- setNames(pal, c("Ancestral", "Derived", "Not extreme"))
	}

	panel <- . %>% 
		plot_ly(x = ~physicalPos, y = ~value) %>%
		add_trace(type = 'scatter',
				mode = 'markers',
				color = ~get(column_color),
				text = ~paste('<b>', rsid,'</b>:', 
							round(value,2),
							'<br><b>', get(column_color), '</b>'),
				hoverinfo = 'text',
				colors = pal,
				showlegend = F) %>%
		add_annotations(
			text = ~unique(POP),
			x = 0.5,
			y = 1,
			yref = "paper",
			xref = "paper",
			yanchor = "bottom",
			showarrow = FALSE,
			font = list(size = 15)
		) %>%
		layout(
			showlegend = FALSE,
			shapes = list(
				list(
					type = "rect",
					x0 = 0,
					x1 = 1,
					xref = "paper",
					y0 = 0, 
					y1 = 16,
					yanchor = 1,
					yref = "paper",
					ysizemode = "pixel",
					fillcolor = toRGB("gray80"),
					line = list(color = "transparent")),
				list(
					type = "line", 
					x0 = 0, 
					x1 = 1, 
					xref = "paper",
					y0 = ~unique(thres), 
					y1 = ~unique(thres), 
					line = list(color = "orange"),
					width = 4, 
					dash = 'dot',
					opacity = 0.4
				),
				list(
					type = "line", 
					x0 = 0, 
					x1 = 1, 
					xref = "paper",
					y0 = ~-unique(thres), 
					y1 = ~-unique(thres), 
					line = list(color = "orange"),
					width = 4, 
					dash = 'dot',
					opacity = 0.4
				)
			),
			xaxis = list(title = 'Position (bp)'),
			yaxis = list(title = varStat)
		)
	
	ByPop_fig <- df %>%
		group_by(POP) %>%
		do(p = panel(.)) %>%
		subplot(
            nrows        = NROW(.),
            shareX       = TRUE,
            shareY       = FALSE,
            which_layout = 1,
            margin       = 0.005,
            heights      = returnHeighPop(popVector)
	)
	return(ByPop_fig)
}

filterThreshold <- function(df_stat_melted,threshold_df,extreme){

	df_stat_melted$thres = 0
	df_stat_melted$thres_ad = df_stat_melted$AD
	out = list()
	for(m in names(threshold_df)){
		tmp = df_stat_melted[POP == m]
		tmp$thres = threshold_df[[m]]
		if(extreme){
			tmp[abs(value) <= threshold_df[[m]]]$thres_ad = 'Not extreme'
			out[[paste0(m)]] = tmp
		}else{
			out[[paste0(m)]] = tmp
		}
	}
	return(rbindlist(out))
}

# To generate the info boxes in summary report
getInfoBoxData <- function(stat,dfVariants){        
	# stat options: c('isafe','ihs','nsl','gwas','dgn', snpeff, clinvar,rgd)
	# Prioretize stats
	keyStats<- data.frame(
		'stat'=c('isafe','ihs','nsl','gwas','snpeff', 'rgd', 'clinvar','dgn'),
		'colNames'=c("isafe","ihsABS","nslABS","NumGWASStudies","rankSnpEFF","RDBranking","ClinVarRank","NofPmids"),'realColNames'=c("isafe","ihs","nsl","NumGWASStudies","effect","RDBranking","ClinicalSignificance","NofPmids"),'order'=c(-1L,-1L,-1L,-1L,1L,1L,1L,-1L))      
	
	orderStats   = append(as.vector(keyStats$colNames[keyStats$stat==stat]),as.vector(keyStats$colNames[keyStats$stat!=stat]))
	direction    = as.vector(keyStats$order[match(orderStats,keyStats$colNames)])
	finalColName = as.vector(keyStats$realColNames[keyStats$stat==stat])

	# Sort data according to stat
	setorderv(dfVariants, orderStats, order = direction,na.last=TRUE)

	# Search max/min value
  topValue  = dfVariants[[orderStats[1]]][1]
  topValues = dfVariants[dfVariants[[orderStats[1]]]==topValue,]
	if(topValue!=0 & !is.na(topValue) & nrow(topValues)>0){
		if(stat=='rgd' & nrow(topValues)>1){
			return(list(topValues$rsid[1],topValues[[finalColName]][1], nrow(topValues),
				paste0(as.vector(topValues$rsid),collapse="%0D%0A")))
		} else {
			return(list(topValues$rsid[1],topValues[[finalColName]][1], nrow(topValues)))
		}
	}
	else{return(NULL)}
}


##################### VARIALBES
# Data for populations
metapopsData <- data.table(
	"meta"=c('AFR', 'EAS' , 'EUR', 'SAS'),
	"colors"=c("#F7F14A", "#33B033", "#5691C4", "#A965BA"))

popsData <- data.table(
	"pops" =c("YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW", 
		"CEU", "TSI", "FIN", "GBR", "IBS",
		"CHB", "JPT", "CHS", "CDX", "KHV", 
		"GIH", "PJL", "BEB", "STU", "ITU"),
	"meta"=c(rep("AFR",7),rep("EUR",5),rep("EAS",5),rep("SAS",5)),
	"colors" =c("#F5F39F","#F0F56B","#F7F14A","#F5F287","#F5EC05",
		"#FAFAC8","#FAFABE","#2D74B2","#ACD1F2","#6FA6D6","#5691C4",
		"#8AB9E3","#33B033","#6DD66D","#37B337","#008F00","#90DE90",
		"#A965BA","#D5A3E3","#8E3EA3","#B79BBD","#C587D6"))


# statCutoff
con <- dbConnect(RMySQL::MySQL(), user='shiny',password='***',
					 dbname='functional',host='localhost')
myQuery <- paste0("SELECT * FROM statCutoff")
statCutoff <- dbGetQuery(con, myQuery) %>% as.data.table
dbDisconnect(con)

# Threshcold by test
ihs_thres = statCutoff[test=='ihs'] %>% select(-c('test'))
nsl_thres = statCutoff[test=='nsl'] %>% select(-c('test'))
isafe_thres = statCutoff[test=='isafe'] %>% select(-c('test'))

# PopsbyMetapopts (KEEP!)
AFRpops<-c('ACB','ASW','ESN','GWD','LWK','MSL','YRI')
EASpops<-c('CDX','CHB','CHS','JPT','KHV')
EURpops<-c('CEU','FIN','GBR','IBS','TSI')
SASpops<-c('BEB','GIH','ITU','PJL','STU')


# VEP
VEPinfo <- data.table(
  "inLabels"= c("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
			  "stop_lost", "start_lost", "bidirectional_gene_fusion", "disruptive_inframe_insertion",
			  "conservative_inframe_insertion", "disruptive_inframe_deletion", "conservative_inframe_deletion",
			  "missense_variant", "splice_region_variant", "start_retained_variant", "stop_retained_variant",
			  "synonymous_variant", "5_prime_UTR_premature_start_codon_gain_variant", "initiator_codon_variant",
			  "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant", "non_coding_transcript_exon_variant",
			  "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant",
			  "intragenic_variant", "intergenic_region"),
  "outLabels"=c("Splice acceptor variant", "Splice donor variant", "Stop gained", "Frameshift variant", "Stop lost",
			  "Start lost", "Bidirectional gene fusion", "Disruptive inframe insertion", "Conservative inframe insertion",
			  "Disruptive inframe deletion", "Conservative inframe deletion", "Missense variant", "Splice region variant",
			  "Start retained variant", "Stop retained variant", "Synonymous variant", "5-UTR premature start codon gain variant",
			  "Initiator codon variant", "5-UTR variant", "3-UTR variant", "Intron variant", "Non-coding transcript exon variant",
			  "Non-coding transcript variant", "Upstream gene variant ", "Downstream gene variant", "Intragenic variant",
			  "Intergenic region"), 'rank'=seq(1,27))

# ClinVar
ClinVarLabels <- data.frame(
	"ClinicalSign"= c('Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic',
	'risk factor','drug response, risk factor','drug response', 'affects','association','protective','Likely benign', 
	'Benign/Likely benign','Benign','other',"Uncertain significance","Conflicting interpretations of pathogenicity",
	"Conflicting interpretations from submitters","no interpretation for the single variant",'association not found',"not provided"), 
		'colors'= c('#D62728','#EF553B', '#E45756','#AB63FA','#625377','#316395','#3B7887', '#3B8775', 
						'#298538','#B6E880','#00CC96','#20701D','#A5C2B2',"#938D8A","#FFA15A", 
						"#FFA15A",'#938D8A','#938D8A','#938D8A'))


#RegulomeDB
regulomeLabels <- c('1a', '1b', '1c', '1d', '1e', '1f', '2a',
				'2b', '2c', '3a', '3b', '4','5','6','7')

ReguColors <- c("#D73027", "#F46D43", "#FDAE61", "#FEE090", "#E0F3F8", 
	"#ABD9E9", "#74ADD1", "#4575B4", "#B3E2CD", "#FDCDAC", "#CBD5E8", 
	"#F4CAE4", "#E6F5C9", "#FFF2AE")


ClinVarPal <- c('#D62728','#EF553B', '#E45756',
		'#AB63FA','#625377','#316395',
		'#3B7887', '#3B8775', '#298538',
		'#B6E880','#00CC96','#20701D','#A5C2B2',
		"#938D8A","#FFA15A", "#FFA15A",'#938D8A','#938D8A','#938D8A')

names(ClinVarPal) <- c('Pathogenic', 
			 'Pathogenic/Likely pathogenic', 
			 'Likely pathogenic',
			 'risk factor',
			 'drug response, risk factor',
			 'drug response', 
			 'affects',
			 'association',
			 'protective',
			 'Likely benign', 
			 'Benign/Likely benign',
			 'Benign',
			 'other',
			 "Uncertain significance",
			 "Conflicting interpretations of pathogenicity",
			 "Conflicting interpretations from submitters",
			 "no interpretation for the single variant",
			 'association not found',
			 "not provided")

# Clocks-GEVA
GEVAclocks <- data.frame(
	"modelClocks" = c("Jnt","Mut","Rec"),
	"plotTitle" = c("Joint clock", "Mutation clock", "Recombination clock"))

# Info Help
help <- read.csv2('helpSteps.csv')

