# print.power.cont.tab
# Print method for objects of class power.cont.tab
# Randy Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created December 22, 2005

print.power.cont.tab <- function(x, ...)
{
  out <- paste("\n\nContingency Table Power Caclulation\n\n",
               "          n1 = ", round(x$n1, digits=1), "\n",
               "          n2 = ", round(x$n2, digits=1), "\n",
               "          p1 = ", round(x$p1, digits=4), "\n",
               "          p2 = ", round(x$p2, digits=4), "\n",
               "          OR = ", round(x$OR, digits=4), "\n",
               "   sig.level = ", round(x$sig.level, digits=4), "\n",
               "       power = ", round(x$power, digits=4), "\n\n\n", sep='')
  cat(out)
}
