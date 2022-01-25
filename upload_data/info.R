suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(DBI))
suppressWarnings(library(RMySQL))

# Connect & get SQL informaiton
sqlToDf <- function(chr, coordStart, coordEnd, dbName){
    con <- dbConnect(RMySQL::MySQL(),
	    user='shiny',
		password='shinyServer?21',
        dbname=dbName,
        host='localhost'
	)
    myQuery <- paste0("SELECT * FROM ", 
    	"chr", chr, 
        " WHERE physicalPos BETWEEN ", coordStart, 
        " AND ", coordEnd
	)
    df <- dbGetQuery(con, myQuery)
    dbDisconnect(con)

    return(df %>% as.data.table)
}

args       = commandArgs(TRUE)
nchr       = as.numeric(gsub("chr","",args[1]))
start      = as.numeric(args[2])
end        = as.numeric(args[3])

#nchr  = 2; start = 135540546 ;end   = 137513561
rsid = fread(paste0(args[4],'rsid.txt'))
names(rsid) = c('physicalPos','rsid')

selscanFiles <- unlist(list.files(path = args[4], pattern=".100bins.norm"))

ihs        = fread(paste0(args[4],selscanFiles[1]))
names(ihs) = c('rsid','physicalPos','freq','ihh1','ihh0','unstandarizediHS','standarizediHS','extremevalue')
ihs = merge(ihs[,c('physicalPos','standarizediHS')],rsid) %>% select(c(physicalPos,rsid,standarizediHS))

nsl        = fread(paste0(args[4],selscanFiles[2]))
names(nsl) = c('rsid','physicalPos','freq','sl1','sl0','unstandarizednSL','standarizednSL','extremevalue')
nsl = merge(nsl[,c('physicalPos','standarizednSL')],rsid) %>% select(c(physicalPos,rsid,standarizednSL))


isafe      = fread(paste0(args[4],"out.iSAFE.out"))[,1:2]
names(isafe)[1] = 'physicalPos'
isafe = merge(isafe,rsid) %>% select(c(physicalPos,rsid,iSAFE))

functional = sqlToDf(local(nchr),local(start),local(end),'functional')
geva = sqlToDf(local(nchr),local(start),local(end),'geva')

df = Reduce(function(x, y) merge(x, y, 
                all=TRUE, by=c("physicalPos", "rsid"), 
                suffixes = c("_x","_y")), 
        list(ihs, nsl, isafe, functional, geva))

df$chr = nchr

fwrite(df,paste0(args[4],"info.txt.gz"))
