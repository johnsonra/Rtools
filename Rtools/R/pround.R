# pround.R
# Randy Johnson
# Laboratory of Genomic Diversity at NCI-Frederick
# SAIC Frederick, Inc
# Last Modified November 20, 2007


pround <- function(p, digits=3, LaTeX=FALSE, full.string=FALSE, equal.length = FALSE)
{
  p <- round(p, digits = digits)

  if(equal.length)
  {
    p.rhs <- sapply(strsplit(as.character(p), '.', fixed = TRUE), function(i) i[2])
    p.len <- nchar(p.rhs)
    suffix <- sapply(digits - p.len, function(i) paste(rep('0', i), collapse = ''))
  }else{
    suffix <- ''
  }
     
  if(LaTeX)
    prefix <- '$<$ 0.'
  else
    prefix <- '< 0.'

  if(full.string)
    return(ifelse(p > 0, paste('p = ', p, suffix, sep=''),
                  paste('p ', prefix, paste(rep(0, digits-1), collapse=''), '1', sep='')))
  else
    return(ifelse(p > 0, paste(p, suffix, sep = ''),
                  paste(prefix, paste(rep(0, digits-1), collapse=''), '1', sep='')))
}
