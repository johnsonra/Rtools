\name{dup.check}
\alias{dup.check}

\title{ Check for duplicate genotypes }
\description{ Check LGD formatted data for duplicate genotyps }

\usage{dup.check(dat, get.list = FALSE)}

\arguments{
  \item{dat}{ LGD formatted data }
  \item{get.list}{ Use this option to have the function return the dup key with hgal numbers }
}

\value{
  \item{num }{Number of duplicate genotypes, (only when get.list==FALSE)}
  \item{dup }{Duplicate key, (only when get.list==TRUE)}
  \item{hgal }{hgal numbers corresponding to the dup key, (only when get.list==TRUE)}
}

\author{ Randy Johnson }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(dat, get.list=FALSE) #dat should be the data in standard LGD format for the gene/SNP in question
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
}
\keyword{ utilities }% at least one, from doc/KEYWORDS
\keyword{ misc }% __ONLY ONE__ keyword per line
