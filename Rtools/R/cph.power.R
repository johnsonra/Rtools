# cph.power.R
# calculate power for a cox proportional hazards model
# Randall Johnson with assistance from Leslie Chinn
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created September 4, 2008
# Last Modified October 29, 2008

cph.power <- function(n1 = NULL, n2 = NULL, n = NULL, prop = NULL,
                      HR = NULL, sig.level = NULL, power = NULL,
                      alternative = "two.sided")
{
  ### checks ###
  # get sample sizes worked out
  if(is.null(n) & !is.null(n1) & !is.null(n2))
    n <- n1 + n2
  if(is.null(n) & !is.null(n2) & !is.null(prop))
    n <- n2 + prop*n2
  if(is.null(n) & !is.null(n1) & !is.null(prop))
    n <- n1/prop + n1
  if(!is.null(n) & !is.null(prop))
    n2 <- n / (1 + prop)
  if(is.null(n1) & !is.null(n) & !is.null(prop))
    n1 <- n / (1/prop + 1)
  if(is.null(prop) & !is.null(n1) & !is.null(n2))
    prop <- n1 / n2

  # if after all this, prop is still NULL, assume prop == 1
  if(is.null(prop))
  {
    NOTE <- "assuming prop == 1"
    prop <- 1
  }

  # check that exactly one of the right things is missing
  if(sum(sapply(list(n, HR, sig.level, power), is.null)) != 1)
    stop("See ?cph.power for details on leaving the correct arguments NULL")

  # number to divide sig.level by in p.body
  denom <- switch(alternative, one.sided = 1, two.sided = 2)

  ### power calculation ###
  p.body <- quote(pnorm((sqrt(prop*n) * abs(HR - 1)) / (prop*HR + 1) -
                        qnorm(1 - sig.level / denom)))

  ### find missing variable ###
  if(is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(1, 1e+07))$root
  else if(is.null(HR))
    HR <- list(uniroot(function(HR) eval(p.body) - power, c(0, 1))$root,
               uniroot(function(HR) eval(p.body) - power, c(1, 1e+04))$root)
  else if(is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(0, 1))$root
  else if(is.null(power))
    power <- eval(p.body)

  # now that I have all of the importaint stuff, make sure I have n1 and n2
  if(is.null(n2))
    n2 <- n / (1 + prop)
  if(is.null(n1))
    n1 <- n / (1/prop + 1)

  # add these to be consistent with power.htest objects
  if(!exists("NOTE"))
    NOTE <- ''
  METHOD <- "Freedman's Method (1982)"

  # return results
  structure(list(n = n, n1 = n1, n2 = n2, prop = prop, HR = HR,
                 sig.level = sig.level, power = power, note = NOTE,
                 method = METHOD), class = "power.htest")
}
