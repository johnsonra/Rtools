\name{smooth.cloud}
\alias{smooth.cloud}

\title{Smooth a cloud of data points}

\description{Returns a vector of points defining a circle around the qth
quantile of the cloud of points}

\usage{smooth.cloud(x, y, q, sf.arc = 1, sf.theta = 1, n = 50, na.rm = FALSE)}

\arguments{
  \item{x}{ Points on the x axis }
  \item{y}{ Points on the y axis }
  \item{q}{ Quanitle of the data to estimate }
  \item{sf.arc}{ Smoothing factor for change in arc length }
  \item{sf.theta}{ Smoothing factor for change in theta }
  \item{n}{ Number of points around the cloud to estimate }
  \item{na.rm}{ Remove NA's from x and y }
}

\details{ x and y are centered at (0,0) and converted to polar
  coordinates. The quantiles are caclulated and back-transformed to
  cartesian coordinates. Weights are normally distributed, but other
  distributions will eventually be included.}

\value{ A data frame with values for the x and y axes }

\author{ Randy Johnson }

\seealso{ wtd.quantile }

\examples{
x <- rnorm(100)
y <- rnorm(100, mean=x)

plot(x, y)
lines(smooth.cloud(x, y, .5), lwd=2)
lines(smooth.cloud(x, y, .95), lwd=2)
}

\keyword{ multivariate }
\keyword{ smooth }
