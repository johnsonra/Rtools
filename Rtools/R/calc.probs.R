# calc.probs.R
# called from polr.power and n.calc
# Calculate probabilities for the treatment group given a cumulative OR and the probabilities for the control group
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 9, 2006

calc.probs <- function(p, OR)
{
  cp <- cumsum(p) # Cumulative Probabilities

  co <- cp / (1-cp) # Cumulative Odds

  coa <- co * OR # Cumulative Odds for Alternate group

  cpa <- coa / (1+coa) # Cumulative Probabilities for Alternate group

  cpa[length(cpa)] <- 1

  return(cpa - c(0, cpa[1:(length(cpa)-1)]))
}
