\name{power.cont.tab}

\alias{power.cont.tab}

\title{Caclulate power for contingency table analyses}
\description{Calculate power, sample size, etc... for contingency table analyses.}

\usage{power.cont.tab(n1 = NULL, n2 = NULL, p1 = NULL, p2 = NULL, sig.level = 0.05, power = NULL, OR = NULL, n2.prop = NULL, alternative = 'two.sided')}

\arguments{
  \item{n1}{Number of cases}
  \item{n2}{Number of controls}
  \item{p1}{Proportion of exposed among cases}
  \item{p2}{Proportion of exposed among controls}
  \item{sig.level}{Significance level, alpha}
  \item{power}{Power, 1 - beta}
  \item{OR}{Odds Ratio}
  \item{n2.prop}{Ratio of controls to cases}
  \item{alternative}{Either 'two.sided' or 'one.sided'}
}

\details{Exactly one of n1, n2, p1, p2, sig.level, and power must be
  NULL. The odds ratio, OR, may be used as a substitute for one, but not both, of the
  proportions. The ratio of controls to cases, n2.prop, may be used as a
  substitute for n1 and n2, or n2, but not n1 only.}

\value{A list of class power.cont.tab with entries corresponding to the input arguements}

\author{Randy Johnson}

\seealso{\code{\link{power.prop.test}}}

\examples{
power.cont.tab(100, 200, .1, .15)
}

\keyword{htest}
