# print.polr.power.R
# Print method for objects of class polr.power
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created March 9, 2006

print.polr.power <- function(x, ...)
{
  out <- paste("\n\nPower of Proportional Odds Logistic Regression\n\n",
               "           n = ", round(x$n, digits=1), "\n",
               "          pc = ", paste(round(x$pc, digits=2), collapse='\t'), "\n",
               "          pt = ", paste(round(x$pt, digits=2), collapse='\t'), "\n\n",
               "          OR = ", round(x$OR, digits=4), "\n",
               "   sig.level = ", round(x$sig.level, digits=4), "\n",
               "       power = ", round(x$power, digits=4), "\n\n\n", sep='')
  cat(out)
}
