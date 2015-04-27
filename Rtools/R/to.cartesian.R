# to.cartessian.R
# Switch from polar to cartesian coordinates
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 23, 2006

to.cartesian <- function(theta, r)
{
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(data.frame(x=x, y=y))
}
