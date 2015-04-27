# MHadd.power.R
# Power calculations for additive genetic model using Mantel-Haenszel's
# chi-squared test for trend
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created February 19, 2010
# Last Modified March 10, 2010

MHadd.power <- function(n1 = NULL, n2 = NULL, f1 = NULL, f2 = NULL, sig.level = 0.05,
                        power = NULL, OR = NULL, n2.prop = NULL, OR.gt.1 = TRUE)
{
  ### Proportions and OR Accounting ###
  # all three should not be specified
  if(!is.null(f2) & !is.null(f1) & !is.null(OR))
  {
    OR <- NULL
    warning("f1, f2, and OR are all specified... OR will be ignored")
  }

  # take care of missing f1 or f2 in presence of OR
  # see calculations in lab book on February 19, 2010 in GWAS section
    if(!is.null(OR) & !is.null(f1)){
      f2 <- f1 / (OR - f1*OR + f1) #f1*OR / (1 - f1 + f1*OR) 
    }else if(!is.null(OR) & !is.null(f2)){
      f1 <- f2*OR / (1 - f2 + f2*OR) #f2 / (OR - f2*OR + f2)
    }

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
  if(sum(sapply(list(n1, n2.prop, f1, f2, power, sig.level), is.null)) != 1)
     stop("exactly one of 'n1', 'n2', 'f1', 'f2', 'power', and 'sig.level' must be NULL")

### calculation of power ###
  p.calc <- quote({
    # sample size
    n <- n1 + n1*n2.prop

    # fill in cells of table
    cell.11 <- (1 - f2)^2 * n1 * n2.prop
    cell.21 <- 2*f2*(1-f2) * n1 * n2.prop
    cell.31 <- f2^2 * n1 * n2.prop
    cell.12 <- (1 - f1)^2 * n1
    cell.22 <- 2*f1*(1-f1) * n1
    cell.32 <- f1^2 * n1

    # calculate allele and disease frequencies
    allele <- (cell.21 + cell.22 + 2*cell.31 + 2*cell.32)/ n
    disease <- (cell.12 + cell.22 + cell.32) / n

    # calculate r
    r.num <- (cell.11*(0-allele)*(0-disease) + cell.12*(0-allele)*(1-disease) +
              cell.21*(1-allele)*(0-disease) + cell.22*(1-allele)*(1-disease) +
              cell.31*(2-allele)*(0-disease) + cell.32*(2-allele)*(1-disease)) / (n - 1)
    r.den <- sqrt((cell.11 + cell.12)*(0-allele)^2 +
                  (cell.21 + cell.22)*(1-allele)^2 +
                  (cell.31 + cell.32)*(2-allele)^2) *
             sqrt((cell.11 + cell.21 + cell.31)*(0-disease)^2 +
                  (cell.12 + cell.22 + cell.32)*(1-disease)^2) / (n - 1)

    r <- r.num / r.den
    
    # find r.alpha
    chisq.alpha <- qchisq(sig.level, df = 1, lower.tail = FALSE)
    r.alpha <- sqrt(chisq.alpha / (n - 1)) # positive

    # transform r.alpha and r
    zp.alpha <- z.prime(r.alpha)
    zp <- z.prime(r)
    zp.se <- 1 / sqrt(n - 3)
    
    # get power :)
    pnorm(-zp.alpha, mean = zp, sd = zp.se) +
      pnorm(zp.alpha, mean = zp, sd = zp.se, lower.tail = FALSE)
  })

### get missing element ###
  if(is.null(n1))    
    n1 <- uniroot(function(n1) eval(p.calc) - power, c(4, 1e+7))$root
  else if(is.null(n2.prop))
    n2.prop <- uniroot(function(n2.prop) eval(p.calc) - power, c(.01, 100))$root
  else if(is.null(f1) & OR.gt.1)
    f1 <- uniroot(function(f1) eval(p.calc) - power, c(f2, 1))$root
  else if(is.null(f2) & OR.gt.1)
    f2 <- uniroot(function(f2) eval(p.calc) - power, c(0, f1))$root
  else if(is.null(f1) & !OR.gt.1)
    f1 <- uniroot(function(f1) eval(p.calc) - power, c(0, f2))$root
  else if(is.null(f2) & !OR.gt.1)
    f2 <- uniroot(function(f2) eval(p.calc) - power, c(f1, 1))$root
  else if(is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.calc) - power, c(0, 1))$root
  else if(is.null(power))
    power <- eval(p.calc)
  else
    stop("internal error")

### Wrap up ###
  if(is.null(OR))
    OR <- ((1 - f2)^2 * 2 * f1 * (1 - f1)) / ((1 - f1)^2 * 2 * f2 * (1 - f2))

  if(is.null(n2))
    n2 <- n1*n2.prop

  return(structure(list(n1 = n1,
                        n2 = n2,
                        f1 = f1,
                        f2 = f2,
                        OR = OR,
                        sig.level = sig.level,
                        power = power),
                   class = 'MHadd.power'))
}

z.prime <- function(r)
  .5*log((1 + r)/(1 - r))
