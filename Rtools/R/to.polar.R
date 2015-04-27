# to.polar.R
# Switch from cartesian to polar coordinates
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 23, 2006

to.polar <- function(x, y)
{
  r <- sqrt(x^2 + y^2)
  theta <- atan2(y, x)
  return(data.frame(theta=theta, r=r))
}
