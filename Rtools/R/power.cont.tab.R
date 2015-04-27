# power.cont.tab.R
# calculate power for contingency tables
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created December 22, 2005
# Last Modified October 29, 2008

power.cont.tab <- function(n1=NULL, n2=NULL, p1=NULL, p2=NULL,
                           sig.level=0.05, power=NULL, OR=NULL,
                           n2.prop=NULL, alternative = 'two.sided')
{
### Proportions and OR Accounting ###
  # all three should not be specified
  if(!is.null(p2) & !is.null(p1) & !is.null(OR))
  {
    OR <- NULL
    warning("p1, p2, and OR are all specified... OR will be ignored")
  }

  # take care of missing p1 or p2 in presence of OR
    if(!is.null(OR) & !is.null(p1))
      p2 <- p.disease(p1, OR)
    else if(!is.null(OR) & !is.null(p2))
      p1 <- p.disease(p2, OR^(-1))

### n1, n2, and n2.prop Accounting ###
  if(is.null(n1))
  {
    if(is.null(n2.prop))
      n2.prop <- 1
    
    if(!is.null(n2))
    {
      warning("n1 is NULL but n2 is specified, ignoring n2")
      n2 <- NULL
    }
  }else if(!is.null(n2))
  {
    if(!is.null(n2.prop))
      warning("n1 and n2 specified, ignoring n2.prop")
              
    n2.prop <- n2/n1
  }  

### Other Accounting ###
  # check that we are missing the correct things
  if(sum(sapply(list(n1, n2.prop, p1, p2, power, sig.level), is.null)) != 1)
     stop("exactly one of 'n1', 'n2', 'p1', 'p2', 'power', and 'sig.level' must be NULL,
           or 'n1' and 'n2' may be NULL iff nequal==TRUE")

### calculation of power ###
  if(alternative == 'two.sided')
  {
    p.calc <- quote({
      cell.11 <- p2 * n1 * n2.prop
      cell.12 <- p1 * n1
      cell.21 <- (1 - p2) * n1 * n2.prop
      cell.22 <- (1 - p1) * n1

      ln.or <- log(p1 * (1-p2) / (p2 * (1-p1)))
      if(abs(ln.or) == Inf)
        return(1)
    
      se.ln.or <- sqrt(1/cell.11 + 1/cell.12 + 1/cell.21 + 1/cell.22)

      pnorm(qnorm(sig.level/2) * se.ln.or, mean=ln.or, sd=se.ln.or) +
        pnorm(qnorm(1 - sig.level/2) * se.ln.or, mean=ln.or, sd=se.ln.or,
              lower.tail=FALSE)
    })
  }else{
    p.calc <- quote({
      cell.11 <- p2 * n1 * n2.prop
      cell.12 <- p1 * n1
      cell.21 <- (1 - p2) * n1 * n2.prop
      cell.22 <- (1 - p1) * n1

      ln.or <- log(p1 * (1-p2) / (p2 * (1-p1)))
      if(abs(ln.or) == Inf)
        return(1)
    
      se.ln.or <- sqrt(1/cell.11 + 1/cell.12 + 1/cell.21 + 1/cell.22)

      if(se.ln.or <= 0)
      {
        eval.sig.level <- 1 - sig.level
        eval.lower.tail <- FALSE
      }else{
        eval.sig.level <- sig.level
        eval.lower.tail <- TRUE
      }
      pnorm(qnorm(eval.sig.level) * se.ln.or, mean=ln.or, sd=se.ln.or,
            lower.tail = eval.lower.tail)
    })
  }

### get missing element ###
  if(is.null(n1))
    n1 <- uniroot(function(n1) eval(p.calc) - power, c(1, 1e+07))$root
  else if(is.null(n2.prop))
    n2.prop <- uniroot(function(n2.prop) eval(p.calc) - power, c(0, 1e+07))$root
  else if(is.null(p1))
    p1 <- uniroot(function(p1) eval(p.calc) - power, c(p2, 1))$root
  else if(is.null(p2))
    p2 <- uniroot(function(p2) eval(p.calc) - power, c(0, p1))$root
  else if(is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.calc) - power, c(0, 1))$root
  else if(is.null(power))
    power <- eval(p.calc)
  else
    stop("internal error")

### Wrap up ###
  if(is.null(OR))
    OR <- p1*(1-p2) / (p2*(1-p1))

  if(is.null(n2))
    n2 <- n1*n2.prop

### Return object ###
  return(structure(list(n1 = n1,
                        n2 = n2,
                        p1 = p1,
                        p2 = p2,
                        OR = OR,
                        sig.level = sig.level,
                        power = power),
                   class="power.cont.tab"))
}
