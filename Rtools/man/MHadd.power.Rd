\name{MHadd.power}

\alias{MHadd.power}

\title{Additive Genetic Model Power using Mantel-Haenszel Test for Trend}

\description{This calculates the power of the Mantel-Haenszel Chi-squared test for trend with the specific application to modeling additively associated SNPs with disease.}

\usage{
MHadd.power(n1 = NULL, n2 = NULL, f1 = NULL, f2 = NULL, sig.level = 0.05, power = NULL, OR = NULL, n2.prop = NULL, OR.gt.1 = TRUE)
}

\arguments{
  \item{n1}{Number of cases}
  \item{n2}{Number of controls}
  \item{f1}{Allele frequency among cases}
  \item{f2}{Allele frequency among controls}
  \item{sig.level}{Significance level, alpha}
  \item{power}{Power, 1 - beta}
  \item{OR}{Odds Ratio -- (change per allele)}
  \item{n2.prop}{Ratio of controls to cases}
  \item{OR.gt.1}{Logical variable indicating wether OR > 1 or OR < 1. This is only used when OR and either f1 or f2 is NULL}
}

\details{Exactly one of n1, n2, f1, f2, sig.level, and power must be
  NULL. The odds ratio, OR, may be used as a substitute for one, but not both, of the
  proportions. The ratio of controls to cases, n2.prop, may be used as a
  substitute for n1 and n2, or n2, but not n1 only.}

\value{A list of class MHadd.power with entries corresponding to the input arguments.}

\references{Need to get some...}

\author{Randy Johnson}

\examples{
# These should both result in roughly the same power
MHadd.power(n1 = 100, n2 = 200, f1 = 0.3, OR = 1.5)

MHadd.power(n1 = 100, f2 = 0.259, OR = 1.5, n2.prop = 2)
}

\keyword{htest}
