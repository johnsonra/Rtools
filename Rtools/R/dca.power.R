# dca.power.R
# Calculate power for Jim's DCA -- Damn Categorical Analysis?
#                                  ...no, Discrete Categorical Analysis
# Randall Johnson
# Created May 20, 2009
# Last Modified July 31, 2009

# need to add some functionality
dca.power <- function(n = NULL, p = NULL, sig.level = 0.05, power = NULL,
                      OR = NULL, n.prop = NULL, n.cat = NULL, n.sim = 1000,
                      codom = FALSE)
{
### n and n.cat Accounting ###
  # number of categories
  n.cat <- length(n)
  
### p and OR Accounting ###
  if(!is.null(OR))
  {
    if(n.cat == 2){
      p <- c(p, p.disease(p, OR))
    }else{
      lnORs <- seq(from = 0, to = log(OR), length = n.cat)
##       lnORs <- 1:n.cat * lnOR # space them appropriately
##       lnORs <- lnORs - mean(lnORs) # center around 0
      ORs <- exp(lnORs)
      p <- p.disease(p, ORs)
    }
  }

### Other Accounting ###
  # very susceptible to user input error at this point... fix it!
  
### Calculation of power ###
  # would be nice to turn off the warnings here...
  do.one.sim <- function(n, p, sig.level)
  {
    # categories and alleles
    cat <- numeric()
    allele <- numeric()

    # assign them for everyone
    for(i in 1:n.cat)
    {
      cat <- c(cat, rep(i-1, n[i]))
      if(!codom){
        allele <- c(allele, rbinom(n[i], 1, p[i]))
      }else{
        allele <- c(allele, rbinom(n[i], 2, p[i]))
      }
    }

    r <- try(cor(cat, allele), silent = TRUE)
    if(class(r) == 'try.error')
      return(FALSE)
    
    p <- pchisq(r * (sum(n) - 1), 1, lower.tail = FALSE)
    return(p <= sig.level)
  }

  # calculate the simulated power
  p.calc <- quote({
    sig <- logical()
    for(i in 1:n.sim)
      sig <- c(sig, do.one.sim(n, p, sig.level))

    sum(sig, na.rm = TRUE) / n.sim
  })

### Get Missing Element ###
    # assumes power is what you want for now
    power <- eval(p.calc)
  
### Wrap up ###

### Return Object ###
  return(structure(list(n = n,
                        p = p,
                        OR = OR,
                        sig.level = sig.level,
                        power = power),
                   class = "dca.power"))
}
