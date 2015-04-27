\name{update.data.frame}

\alias{update.data.frame}

\title{Update a data frame}

\description{Update a data frame similar to PROC UPDATE in SAS}

\usage{update.data.frame(object, formula, ..., newdata)}

\arguments{
  \item{object}{Data frame to be updated.}
  \item{formula}{Two sided formula with the id variable(s) on the left hand side and the variables to be updated on the right hand side (see Details below).}
  \item{\dots}{Included to be consistent with update(), but this is effectively ignored.}
  \item{newdata}{Data frame with updated data.}
}

\details{The formula argument should be specified as in the example section. Each id variable must be in both the original (object) and the new (newdata) data frames, but the variables to be updated need not be in the original data frame. All ids not present in the new data frame will be replaced with missing values. Each id or combination of ids must be unique in the new data frame as well as for every row in the original data frame that is to be updated. Any variable in the original data frame that is not included in the formula will be untouched, and any variable in the new data frame will be ignored.}

\value{A data frame with all data object replaced by corresponding data in newdata. Data that has no replacement remains unchanged. The output data frame is sorted on the id variables.}

\author{Randall Johnson}

\examples{
set.seed(29084)

x1 <- rnorm(30)
x2 <- rnorm(30)

dat.ori <- data.frame(id1 = rep(1:10, each = 3),
                      id2 = rep(c('a', 'b', 'c'), 10),
                      x1 = x1,
                      x2 = x2,
                      y = 2*x1 - 2*x2 + rnorm(30))

x1 <- x1[-c(1:10*3)]
x2 <- x2[-c(1:10*3)]

dat.new <- data.frame(id1 = rep(1:10, each = 2),
                      id2 = rep(c('a', 'b'), 10),
                      x2 = x1 + rnorm(20, sd = 2),
                      y = 2*x1 - 2*x2 + rnorm(20))

dat.updated <- update(dat.ori, id1 + id2 ~ x2 + y, newdata = dat.new)

dat.new$z <- rbinom(20, size = 2, prob = .5)

dat.updated <- update(dat.ori, id1 + id2 ~ z, newdata = dat.new)
}

\keyword{manip}
\keyword{methods}
