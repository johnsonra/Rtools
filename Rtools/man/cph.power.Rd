\name{cph.power}

\alias{cph.power}

\title{Power and Sample Size estimation Under Proportional Hazards Model}

\description{Estimation of power and sample size given the number of events, significance level, frequency of exposure, and hazard ratio.}

\usage{
cph.power(n1 = NULL, n2 = NULL, n = NULL, prop = NULL, HR = NULL, sig.level = NULL, power = NULL, alternative = "two.sided")
}

\arguments{
  \item{n1}{Number of cases}
  \item{n2}{Number of controls}
  \item{n}{Total numer of events}
  \item{prop}{Proportion of cases to controls (ie n1 / n2 or phenotype frequency)}
  \item{HR}{Hazard Ratio}
  \item{sig.level}{Significance level}
  \item{power}{Power}
  \item{alternative}{Alternative hypothesis ("one.sided" or "two.sided")}
}

\details{Exactly one of 'num.events', 'freq', 'HR', 'sig.level', or 'power' must be NULL}

\value{An object of class htest}

\references{Freedman 1982 in Statistics in Medicine, 11, 179-91}

\author{Randall Johnson and Leslie Chinn}

\examples{
cph.power(n = 172, prop = 1, HR = 0.7, sig.level = 0.05)
}

\keyword{htest}
