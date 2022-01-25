# This script is adapted from https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/blob/master/interpolate_maps.R

suppressMessages(suppressWarnings(library(intervals)))
suppressMessages(library(data.table))

getRate = function(map, intI) {
        # get recombination rate in cM/Mb for a given map and interval
        rr1 <- approx( map$pos, map$cM, xout=intI[,1] )
        rr2 <- approx( map$pos, map$cM, xout=intI[,2] )
        rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # cM/Mb

        cm <- approx(map$pos, map$cM, xout=intI[,1]) # cM

        # Store both values in a list
        output <- list()
        output[['rate']] <- rate
        output[['cM']] <- cm
        return(output)
}

args    = commandArgs(TRUE)

pos     = fread(args[1])
pos[,4] = pos[,4] + 1
nchr    = unique(pos[,1])
nchr    = as.numeric(gsub('chr','',nchr))
map     = fread(paste0('/home/pophumanvar/phv_pipeline/rec_map/sexavg_chr',nchr,'.txt'))

tmp     = getRate(map,as.matrix(pos[,3:4]))

df      = cbind(pos[,1],paste0('P',seq(0,nrow(pos)-1)),tmp[[2]][[2]],pos[,3])
fwrite(df,args[2],sep='\t',col.names=F)
