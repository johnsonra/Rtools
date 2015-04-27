# polr.power.R
# Calculate power for Proportional Odds Logistic Regression
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick
# Created March 8, 2006
# Last Modified March 9, 2006

polr.power <- function(n=NULL, pc=NULL, pt=NULL, sig.level=0.05, power=0.8, OR=NULL)
{
### check that pc and pt sum to 1 ###
  if((!is.null(pt) & sum(pt) != 1) | (!is.null(pc) & sum(pc) != 1))
    stop("pc and pt must sum to 1 if specified")
  
### OR / pc / pt Accounting ###
  if(!is.null(pt) & !is.null(OR) & is.null(pc))
    pc <- calc.probs(pt, 1/OR)
  else if(!is.null(pc) & !is.null(pt) & is.null(OR))
  {
    warning("Estimating OR from pc and pt... pt may be altered in output")
    cpc <- cumsum(pc)
    cpt <- cumsum(pt)
    OR <- mean((cpc / (1-cpc)) / (cpt / (1-cpt)), na.rm=TRUE)
  }else if(is.null(pc) & (is.null(OR) | is.null(pt)))
    stop("pc must be specified if either OR or pt is missing")
  else if(!is.null(pc) & !is.null(pt) & !is.null(OR))
    warning("OR, pc and pt need not be specified together... pt may be altered in output")

### Other Accounting ###
  if(sum(sapply(list(n, pc, sig.level, power, OR), is.null)) != 1)
    stop("Exactly one of n, pc, sig.level, power, or OR may be null")
  
### Get answer ###
  if(is.null(n))
    n <- n.calc(pc, OR, sig.level, power)
  else if(is.null(OR))
    OR <- uniroot(function(x){n.calc(pc, OR=x, sig.level, power) - n}, c(1, 1e+07))$root
  else if(is.null(sig.level))
    sig.level <- uniroot(function(x){n.calc(pc, OR, sig.level=x, power) - n}, c(0, 1))$root
  else if(is.null(power))
  {
    power <- try(uniroot(function(x){n - n.calc(pc, OR, sig.level, power=x)}, c(sig.level, 1))$root, silent=TRUE)
    if(class(power) == 'try-error')
      power <- sig.level
  }
  else
    stop("internal error")
  
### Return object ###
  return(structure(list(n = n,
                        pc = pc,
                        pt = calc.probs(pc, OR),
                        OR = OR,
                        sig.level = sig.level,
                        power = power),
                   class="polr.power"))
}
