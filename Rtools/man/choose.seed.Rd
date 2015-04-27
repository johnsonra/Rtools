\name{choose.seed}
\alias{choose.seed}

\title{ Choose a random seed }

\description{This will choose a random seed of a specified magnitude... the implementation is a bit outdated.}

\usage{choose.seed(len = 9)}

\arguments{ \item{len}{ Length of the seed (determining the magnitude of the integer) } }

\value{ A random integer of a given magnitude is returned }

\author{ Randy Johnson }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(len=9)
{
    x = sample(c(0:9), len, replace=TRUE)
    return(sum(x*(rep(10,len)^c(0:(len-1)))))
  }
}

\keyword{ utilities }% at least one, from doc/KEYWORDS
\keyword{ misc }% __ONLY ONE__ keyword per line
