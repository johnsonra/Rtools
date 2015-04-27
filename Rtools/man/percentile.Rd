\name{percentile}

\alias{percentile}

\title{ Calculate a percentile }

\description{ Calculate the percentile of a number given an emperical distribution }

\usage{ percentile(x, distn) }

\arguments{
  \item{x}{ A number of vector of numbers to calculate percentiles for }
  \item{distn}{ A vector with an emperical distribution }
}

\value{ The percentile or a vector of percentiles is returned }

\author{ Randy Johnson }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, distn)
{
    distn = subset(distn, !is.na(distn))
    n = length(distn)
    
    if(n<2 | is.null(distn) | max(distn)==min(distn))
    { 
        warning("Distribution vector is too short, NULL, or has only one value in it.")
        return(NA)
    }
    if(is.null(x))
    { 
        warning("Vector of quantiles is NULL.")
        return(NA)
    }

    y = findInterval(x, distn[order(distn)]) / n
    
    return(y)
  }
}

\keyword{ utilities }% at least one, from doc/KEYWORDS
\keyword{ misc }% __ONLY ONE__ keyword per line
