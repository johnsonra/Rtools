# p.disease.R
# helper function for power.cont.tab
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created December 22, 2005
# Last Modified January 27, 2006

p.disease <- function(p, OR) # see page 56 in lab book starting August 2004
{
  return(OR*p / (1 + OR*p - p))
}
