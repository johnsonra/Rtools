#Randy Johnson
#R code for checking for duplicate genotypes
#Created Sept 7, 2004
#Last Modified Sept 7, 2004

dup.check = function(dat, get.list=FALSE) #dat should be the data in standard LGD format for the gene/SNP in question
{
    dat = unique(dat) # the rest of the dat data part is to remove duplicates with different genotypes (ie one person with two different genotypes)
    dat = dat[order(dat$hgal),]
    dat$dup = rep(0, length(dat$hgal))
    for(i in c(1:(length(dat$hgal)-1)))
    {
        dat$dup[i] = ifelse(!is.na(dat$hgal[i]) & !is.na(dat$hgal[i+1]) & dat$hgal[i]==dat$hgal[i+1], 1, dat$dup[i])
    }
    for(i in c(2:length(dat$hgal)))
    {
        dat$dup[i] = ifelse(!is.na(dat$hgal[i]) & !is.na(dat$hgal[i-1]) & dat$hgal[i]==dat$hgal[i-1], 1, dat$dup[i])
    }
    if(get.list){return(subset(dat, select=c('dup', 'hgal')))
    }else{return(sum(dat$dup))}
}
