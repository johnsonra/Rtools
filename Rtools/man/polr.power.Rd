\name{polr.power}
\alias{polr.power}

\title{Power and Sample Size Calculations for Proportional Odds Logistic
Regression}

\description{Calculates power, sample size, OR's, etc... for
  proportional odds logistic regression.}

\usage{polr.power(n = NULL, pc = NULL, pt = NULL, sig.level = 0.05, power = 0.8, OR = NULL)}

\arguments{
  \item{n}{Sample size per group}
  \item{pc}{Probabilities for control group}
  \item{pt}{Probabilities for treatment group}
  \item{sig.level}{Significance level}
  \item{power}{Power}
  \item{OR}{Cumulative Odds Ratio}
}

\value{
  A polr.power object: List of 6 objects
  \item{n}{Sample size per group}
  \item{pc}{Probabilities for control group}
  \item{pt}{Probabilities for treatment group}
  \item{sig.level}{Significance level}
  \item{power}{Power}
  \item{OR}{Cumulative Odds Ratio}
}

%\references{}

\author{Randy Johnson}

\seealso{\code{polr}}

\examples{
pc <- c(0.10, 0.20, 0.30, 0.25, 0.15)

# how big would my sample have to be to detect
# a cumulative odds ratio of at least 2 with
# 80\% power at a significance level of 0.05?

polr.power(pc=pc, OR=2)

# If I can only get 50 samples per group, what
# is the largest cumulative odds ratio I can
# expect to detect with 80\% power at a significance
# level of 0.01?

polr.power(n=50, pc=pc, sig.level=0.01)
}

\keyword{multivariate}
\keyword{htest}
