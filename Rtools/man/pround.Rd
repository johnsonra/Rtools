\name{pround}

\alias{pround}

\title{Round p-values}

\description{Round p-values with appropriate output for small p-values.}

\usage{pround(p, digits = 3, LaTeX = FALSE, full.string=FALSE, equal.length = FALSE)}

\arguments{
  \item{p}{Vector of p-values}
  \item{digits}{num Number of digits to round to}
  \item{LaTeX}{logical When TRUE '<' is output as '$<$' to compile
    properly in LaTeX}
  \item{full.string}{logical When TRUE 'p = ' or 'p < ' is added to the
    string}
  \item{equal.length}{make all p-values the same length (eg `1.00'
  instead of `1')}
}

\value{Rounded p-values}

\author{Randy Johnson}

\examples{
pround(c(0, 0.54587379, 0.837462, 0.48579384, 0.0000003, 0.0002, 0.2938472))
}

\keyword{manip}
