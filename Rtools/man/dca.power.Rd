\name{dca.power}

\alias{dca.power}

\title{DCA Power}

\description{Calculate the power for Jim's DCA.}

\usage{
dca.power(n = NULL, p = NULL, sig.level = 0.05, power = NULL, OR = NULL, n.prop = NULL, n.cat = NULL, n.sim = 1000, codom = FALSE)
}

\arguments{
  \item{n}{Vector of number of individuals in each category (ordered).}
  \item{p}{Phenotype frequency of the phenotype associated with the disease (eg AG/GG for a dominant model) when codom is FALSE, or the minor allele frequency when codom is TRUE, in the reference group}
  \item{sig.level}{Significance level}
  \item{power}{Power -- ignored for now}
  \item{OR}{Odds ratio for comparison of one group to the next highest group}
  \item{n.prop}{Porportion of individuals in one group to the next -- ignored for now}
  \item{n.cat}{Number of categories -- ignored for now}
  \item{n.sim}{Number of repetitions to run.}
  \item{codom}{Logical, indicating if it is a codominant genetic model (assumed dominant or recessive when FALSE)}
}

\value{
A list of class `dca.power' with the following components:
  \item{n}{Number of individuals in each category}
  \item{p}{Frequency of the minor allele in each category}
  \item{OR}{Odds ratio for comparison of one group to the next highest group}
  \item{sig.level}{Significance level}
  \item{power}{Power of the test}
}

\author{Randall Johnson}

\examples{
# need to add one
}

\keyword{htest}
