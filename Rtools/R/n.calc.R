# n.calc.R
# Called from polr.power
# Calculate sample size per group given control probabilities, cumulative OR, significance level, and power
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 9, 2006

n.calc <- function(pc, OR, sig.level, power)
{
  pt <- calc.probs(pc, OR)
  
  pibar <- sapply(1:length(pc), function(i){mean(c(pc[i], pt[i]))})
    
  Z <- qnorm(1-sig.level/2) + qnorm(power)

  n <- 6 * Z^2 / (log(OR)^2) / (1 - sum(pibar^3))
  if(n == Inf)
    return(1e+07)
  else
    return(n)
}
