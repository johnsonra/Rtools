# print.MHadd.power.R
# Print method for objects of class MHadd.power
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Adapted from print.power.cont.tab February 20, 2010

print.MHadd.power <- function(x, ...)
{
  out <- paste("\n\nAdditive Genetic Model Power via Mantel-Haenszel Trend Test\n\n",
               "          n1 = ", round(x$n1, digits=1), "\n",
               "          n2 = ", round(x$n2, digits=1), "\n",
               "          f1 = ", round(x$f1, digits=4), "\n",
               "          f2 = ", round(x$f2, digits=4), "\n",
               "          OR = ", round(x$OR, digits=4), "\n",
               "   sig.level = ", round(x$sig.level, digits=4), "\n",
               "       power = ", round(x$power, digits=4), "\n\n\n", sep='')
  cat(out)
}
